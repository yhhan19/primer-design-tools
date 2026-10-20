#include "pipeline.hpp"

void PDRStage::run(PipelineContext& ctx) {
    srand(ctx.args.seed);

    const bool midpoint_search = ctx.args.pdr_search == "midpoint";
    std::cout << "PDR search mode: "
              << (midpoint_search
                      ? "midpoint-constrained variable amplicon length"
                      : "fixed amplicon length")
              << "\n";

    RiskOptimizer ro(ctx.labels,
                     ctx.sequences, 
                     ctx.args.len_PDR, 
                     midpoint_search ? ctx.args.len_amp_min
                                     : ctx.args.len_amp,
                     ctx.args.len_amp,
                     ctx.args.len_fu,
                     ctx.args.len_ru,
                     ctx.args);

    auto PDR = ro.search(0, ctx.args.u_max, ctx.args.u_min);

    if (PDR.empty()) {
        throw std::runtime_error(
            "PDR optimization returned no feasible primer-design regions");
    }

    ro.validate_PDR(PDR);
    
    risk_t score = ro.score(PDR, ALPHA);
    std::cout << "loss: " << score << std::endl;

    std::ofstream fout(ctx.args.output_file + "." + short_name());
    fout << display_PDR(PDR) << std::endl;
    fout.close();

    ctx.pdr_regions = PDR;
}

void PrimerSelStage::run(PipelineContext& ctx) {
    std::cout << "Starting primer selection for "
              << ctx.pdr_regions.size() / 2 << " PDRs\n";

    Thal::init(std::string(PRIMER3_PATH) + "/src/primer3_config",
               ctx.args.mv, ctx.args.dv, ctx.args.dntp,
               ctx.args.dna_conc, ctx.args.temp);

    // ==============================================================
    // 1. Primer3 global settings
    // ==============================================================
    p3_global_settings *pa = p3_create_global_settings();

    p3_set_gs_primer_opt_size(pa, ctx.args.p3_opt_size);
    p3_set_gs_primer_min_size(pa, ctx.args.p3_min_size);
    p3_set_gs_primer_max_size(pa, ctx.args.p3_max_size);

    p3_set_gs_primer_opt_tm(pa, ctx.args.p3_opt_tm);
    p3_set_gs_primer_min_tm(pa, ctx.args.p3_min_tm);
    p3_set_gs_primer_max_tm(pa, ctx.args.p3_max_tm);

    p3_set_gs_primer_min_gc(pa, ctx.args.p3_min_gc);
    p3_set_gs_primer_max_gc(pa, ctx.args.p3_max_gc);

    pa->p_args.salt_conc     = ctx.args.mv;
    pa->p_args.divalent_conc = ctx.args.dv;
    pa->p_args.dntp_conc     = ctx.args.dntp;
    pa->p_args.dna_conc      = ctx.args.dna_conc;

    // Primer3 validates SEQUENCE_PRIMER_PAIR_OK_REGION_LIST against the
    // product interval even in pick_primer_list mode. We do not construct
    // pairs here, so a biological product-size interval would incorrectly
    // reject valid independent oligo regions (especially a minimum-length
    // midpoint PDR). The PDR optimizer already enforces the requested
    // amplicon range. Use the broadest legal interval for validation:
    // Primer3 requires its lower bound to be at least PRIMER_MAX_SIZE.
    pa->pr_min[0] = static_cast<int>(ctx.args.p3_max_size);
    pa->pr_max[0] = static_cast<int>(std::min<std::size_t>(
        ctx.tmpl.size(),
        static_cast<std::size_t>(std::numeric_limits<int>::max())));
    pa->num_intervals = 1;

    // Generate complete oligo lists first. Retained candidates are explicitly
    // paired and validated before downstream stages.  Assign the enum directly: the
    // string task setter treats pick_primer_list as a legacy shortcut and
    // may also enable internal oligos.
    pa->primer_task         = pick_primer_list;
    pa->pick_left_primer    = 1;
    pa->pick_right_primer   = 1;
    pa->pick_internal_oligo = 0;

    std::cout << "Primer3 task: pick_primer_list"
              << "  |  Template length: "
              << ctx.tmpl.size()
              << " bp\n";


    // ==============================================================
    // 2. Primer3 sequence arguments
    // ==============================================================
    seq_args *sa = create_seq_arg();

    sa->sequence =
        strdup(ctx.tmpl.c_str());

    sa->sequence_name =
        strdup("my_template");

    sa->ok_regions.any_left  = 0;
    sa->ok_regions.any_right = 0;
    sa->ok_regions.any_pair  = 0;
    sa->ok_regions.count     = 1;

    sa->ok_regions.left_pairs[0][1] =
        ctx.args.len_PDR;

    sa->ok_regions.right_pairs[0][1] =
        ctx.args.len_PDR;


    // ==============================================================
    // 3. Ranking rule for variant robustness
    //
    // Lower is better.
    //
    // Priority:
    //   1. weighted mismatch
    //   2. 3'-end mismatch rate
    //   3. genome mismatch rate
    //   4. ordinary per-base mismatch rate
    // ==============================================================
    auto variant_better =
        [](const Oligo& a, const Oligo& b) {

            if (a.weighted_mismatch !=
                b.weighted_mismatch) {

                return a.weighted_mismatch <
                       b.weighted_mismatch;
            }

            if (a.mismatch_3prime_rate !=
                b.mismatch_3prime_rate) {

                return a.mismatch_3prime_rate <
                       b.mismatch_3prime_rate;
            }

            if (a.genome_mismatch_rate !=
                b.genome_mismatch_rate) {

                return a.genome_mismatch_rate <
                       b.genome_mismatch_rate;
            }

            return a.mismatch_rate <
                   b.mismatch_rate;
        };

    // ==============================================================
    // 6. Run Primer3 for each PDR pair
    // ==============================================================
    struct KeptPrimerSummary {
        int pdr_index;
        index_t left_start;
        index_t right_start;
        std::size_t left_raw;
        std::size_t right_raw;
        std::size_t left_paired = 0;
        std::size_t right_paired = 0;
        PairingStats pairing;
        std::size_t msa_kept = 0;
    };
    std::vector<KeptPrimerSummary> kept_summary;
    kept_summary.reserve(ctx.pdr_regions.size() / 2);

    for (std::size_t i = 0;
         i < ctx.pdr_regions.size();
         i += 2) {

        const int region_idx =
            static_cast<int>(i / 2);


        sa->ok_regions.left_pairs[0][0] =
            ctx.pdr_regions[i];

        // _adjust_seq_args() may optimize these intervals in place, so
        // restore the requested PDR lengths before every Primer3 call.
        sa->ok_regions.left_pairs[0][1] =
            ctx.args.len_PDR;

        sa->ok_regions.right_pairs[0][0] =
            ctx.pdr_regions[i + 1];

        sa->ok_regions.right_pairs[0][1] =
            ctx.args.len_PDR;


        // ----------------------------------------------------------
        // Run Primer3
        // ----------------------------------------------------------
        p3retval *retval =
            choose_primers(pa, sa);


        if (retval == nullptr) {

            std::cerr
                << "ERROR: choose_primers() returned NULL "
                << "for region "
                << region_idx
                << "\n";

            destroy_seq_args(sa);
            p3_destroy_global_settings(pa);

            return;
        }


        if (retval->glob_err.data) {

            std::cerr
                << "Global error: "
                << retval->glob_err.data
                << "\n";
        }


        if (retval->per_sequence_err.data) {

            std::cerr
                << "Sequence error: "
                << retval->per_sequence_err.data
                << "\n";
        }

        if (retval->warnings.data) {
            std::cerr
                << "Primer3 warning: "
                << retval->warnings.data
                << "\n";
        }


        // ==========================================================
        // 7. Extract Primer3 candidates
        // ==========================================================
        PrimerOutput out =
            extract_all(retval, sa);

        out.pdr_left =
            ctx.pdr_regions[i];

        out.pdr_right =
            ctx.pdr_regions[i + 1];

        out.segment_id =
            ctx.current_seg_id();


        const std::size_t left_before =
            out.left_oligos.size();

        const std::size_t right_before =
            out.right_oligos.size();


        // ==========================================================
        // 8. Variant-aware scoring
        // ==========================================================

        // Left candidates
        for (Oligo& o : out.left_oligos) {

            score_primer_variants(
                o,
                ctx.sequences,
                false
            );
        }


        // Right candidates
        for (Oligo& o : out.right_oligos) {

            score_primer_variants(
                o,
                ctx.sequences,
                true
            );
        }


        // ==========================================================
        // 9. Rank candidates for reporting; retain all before pairing
        // ==========================================================
        std::stable_sort(
            out.left_oligos.begin(),
            out.left_oligos.end(),
            variant_better
        );

        std::stable_sort(
            out.right_oligos.begin(),
            out.right_oligos.end(),
            variant_better
        );

        if (out.left_oligos.empty() || out.right_oligos.empty()) {
            destroy_p3retval(retval);
            destroy_seq_args(sa);
            p3_destroy_global_settings(pa);
            throw std::runtime_error(
                "Primer3 returned no candidates on at least one side "
                "for PDR " + std::to_string(region_idx));
        }

        auto report_shortage = [&](const char* side, std::size_t count) {
            if (count < ctx.args.num_return) {
                std::cerr << "PDR " << region_idx << " " << side
                          << " available candidates: " << count
                          << " < requested " << ctx.args.num_return
                          << "\n";
            }
        };
        report_shortage("LEFT", out.left_oligos.size());
        report_shortage("RIGHT", out.right_oligos.size());

        kept_summary.push_back({
            region_idx,
            ctx.pdr_regions[i],
            ctx.pdr_regions[i + 1],
            left_before,
            right_before
        });


        // ==========================================================
        // 11. Save all ranked candidates
        // ==========================================================
        const auto pairing_stats = build_candidate_pairs(out, ctx.tmpl, ctx.args);
        auto& summary = kept_summary.back();
        summary.pairing = pairing_stats;
        if (out.pairs.empty()) {
            destroy_p3retval(retval);
            destroy_seq_args(sa);
            p3_destroy_global_settings(pa);
            throw std::runtime_error("No valid forward/reverse candidate pair for PDR " +
                std::to_string(region_idx) + " in segment " + std::to_string(out.segment_id));
        }
        score_candidate_pairs_msa(out, ctx.sequences, ctx.args);
        const auto ranked_pairs = out.pairs;  // Retain a full audit of accepted and dropped pairs.
        summary.msa_kept = retain_candidate_pairs_msa(out, ctx.args);
        summary.left_paired = out.left_oligos.size();
        summary.right_paired = out.right_oligos.size();
        const std::string msa_path = ctx.args.output_file + ".segment" +
            std::to_string(out.segment_id) + ".pdr" + std::to_string(region_idx) + ".pair-msa.tsv";
        std::ofstream msa_file(msa_path);
        if (!msa_file) throw std::runtime_error("Cannot write " + msa_path);
        msa_file << "rank\tleft_start\tright_start\tleft_reference_sequence\tright_reference_sequence"
                 << "\tjoint_coverage\tjoint_perfect_match\tjoint_3prime_match\tweighted_mismatch\tunknown_fraction\tmsa_sequences\tkept\n";
        msa_file << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (std::size_t rank = 0; rank < ranked_pairs.size(); ++rank) {
            const auto& pair = ranked_pairs[rank]; const auto& m = pair.msa;
            msa_file << rank+1 << '\t' << pair.left.start << '\t' << pair.right.start
                     << '\t' << pair.left.seq << '\t' << pair.right.seq << '\t' << m.joint_coverage
                     << '\t' << m.joint_perfect_match << '\t' << m.joint_3prime_match
                     << '\t' << m.weighted_mismatch << '\t' << m.unknown_fraction
                     << '\t' << m.sequence_count << '\t' << (rank < summary.msa_kept) << '\n';
        }
        msa_file.close();
        if (!msa_file) throw std::runtime_error("Failed writing " + msa_path);
        const std::string pair_path = ctx.args.output_file + ".segment" +
            std::to_string(out.segment_id) + ".pdr" + std::to_string(region_idx) + ".pairs.tsv";
        std::ofstream pair_file(pair_path);
        if (!pair_file) throw std::runtime_error("Cannot write " + pair_path);
        pair_file << "pair_index\tleft_start\tright_start\tleft_sequence\tright_reference_sequence\tproduct_bp\n";
        for (std::size_t pair_id = 0; pair_id < out.pairs.size(); ++pair_id) {
            const auto& pair = out.pairs[pair_id];
            pair_file << pair_id << '\t' << pair.left.start << '\t' << pair.right.start
                      << '\t' << pair.left.seq << '\t' << pair.right.seq << '\t' << pair.product_size << '\n';
        }
        pair_file.close();
        if (!pair_file) throw std::runtime_error("Failed writing " + pair_path);
        ctx.candidate_primers.push_back(std::move(out));


        destroy_p3retval(retval);
    }


    // ==============================================================
    // 13. Final retained-primer summary
    // ==============================================================
    constexpr int SUMMARY_WIDTH = 134;
    std::size_t total_left_raw = 0, total_right_raw = 0;
    std::size_t total_left_paired = 0, total_right_paired = 0;
    std::size_t total_combinations = 0, total_checked = 0;
    std::size_t total_p3_rejected = 0, total_pairs = 0, total_msa_kept = 0;
    std::cout << "\nCANDIDATE PAIR SUMMARY (after Primer3 checks and MSA filtering)\n"
              << "In pairs/raw = distinct primers present in retained MSA-ranked pairs / original candidates.\n"
              << std::string(SUMMARY_WIDTH, '-') << '\n' << std::right
              << std::setw(8) << "PDR" << std::setw(12) << "Left PDR"
              << std::setw(12) << "Right PDR" << std::setw(20) << "Left in pairs/raw"
              << std::setw(20) << "Right in pairs/raw" << std::setw(12) << "Combinations"
              << std::setw(14) << "P3 checked" << std::setw(12) << "P3 rejected"
              << std::setw(12) << "P3 valid" << std::setw(12) << "MSA kept" << '\n'
              << std::string(SUMMARY_WIDTH, '-') << '\n';
    for (const auto& row : kept_summary) {
        const auto p3_rejected = row.pairing.checked - row.pairing.accepted;
        total_left_raw += row.left_raw; total_right_raw += row.right_raw;
        total_left_paired += row.left_paired; total_right_paired += row.right_paired;
        total_combinations += row.pairing.combinations;
        total_checked += row.pairing.checked;
        total_p3_rejected += p3_rejected; total_pairs += row.pairing.accepted;
        total_msa_kept += row.msa_kept;
        std::cout << std::setw(8) << row.pdr_index << std::setw(12) << row.left_start
                  << std::setw(12) << row.right_start
                  << std::setw(20) << (std::to_string(row.left_paired) + "/" + std::to_string(row.left_raw))
                  << std::setw(20) << (std::to_string(row.right_paired) + "/" + std::to_string(row.right_raw))
                  << std::setw(12) << row.pairing.combinations
                  << std::setw(14) << row.pairing.checked << std::setw(12) << p3_rejected
                  << std::setw(12) << row.pairing.accepted << std::setw(12) << row.msa_kept << '\n';
    }
    std::cout << std::string(SUMMARY_WIDTH, '-') << '\n'
              << std::setw(32) << "TOTAL"
              << std::setw(20) << (std::to_string(total_left_paired) + "/" + std::to_string(total_left_raw))
              << std::setw(20) << (std::to_string(total_right_paired) + "/" + std::to_string(total_right_raw))
              << std::setw(12) << total_combinations
              << std::setw(14) << total_checked << std::setw(12) << total_p3_rejected
              << std::setw(12) << total_pairs << std::setw(12) << total_msa_kept << '\n'
              << std::string(SUMMARY_WIDTH, '-') << '\n'
              << "P3 checked: combinations submitted to Primer3; no geometry/product-length prefilter.\n"
              << "P3 rejected: Primer3 returned no acceptable pair.\n"
              << "MSA kept: ranked pairs retained after percentage/minimum filtering, including boundary ties.\n"
              << "Primers without a retained pair are excluded from downstream candidate lists.\n";

    // ==============================================================
    // 14. Cleanup
    // ==============================================================
    destroy_seq_args(sa);
    p3_destroy_global_settings(pa);
}

void OffTargetStage::run(PipelineContext& ctx) {
    if (ctx.args.ref_file == "") {
        std::cout << "No ref file provided, skip off target search\n";
        for (auto& seg : ctx.segments)
            seg.filtered_primers = seg.candidate_primers;
        return;
    }

    Thal::init(std::string(PRIMER3_PATH) + "/src/primer3_config",
               ctx.args.mv, ctx.args.dv, ctx.args.dntp,
               ctx.args.dna_conc, ctx.args.temp);

    // flatten all candidate_primers from all segments — track offsets
    std::vector<PrimerOutput> all_candidates;
    std::vector<std::size_t>  seg_offsets;  // seg_offsets[s] = start index in all_candidates

    for (auto& seg : ctx.segments) {
        seg_offsets.push_back(all_candidates.size());
        for (auto& po : seg.candidate_primers)
            all_candidates.push_back(po);
    }

    std::cout << "Off-target search: " << all_candidates.size()
              << " PDRs across " << ctx.segments.size() << " segments\n";

    // single automaton + single genome scan
    std::vector<std::string> labels, sequences;
    convert(all_candidates, labels, sequences);

    Automaton ac(labels, sequences, ctx.args.kmer_len);
    std::vector<PrimerOutput> all_filtered;
    if (ctx.args.off_target_mode == "single") {
        auto results = ac.search(ctx.args.ref_file, ctx.args.kmer_len - 1,
            ctx.args.threshold, ctx.args.dg_thres, ctx.args.chunk_size,
            ctx.args.block_size, ctx.args.nthreads);
        write_results(ctx.args.output_file + "." + short_name(), results, labels);
        all_filtered = filterByDG_relax(results, ctx.args.dg_thres,
            all_candidates, ctx.args.num_return);
    } else {
        std::vector<AmpliconPrimerPair> definitions;
        for (std::size_t oi = 0; oi < all_candidates.size(); ++oi) {
            const auto& candidates = all_candidates[oi];
            if (candidates.pairs.empty())
                throw std::runtime_error("Paired off-target screening requires retained candidate pairs");
            for (std::size_t pi = 0; pi < candidates.pairs.size(); ++pi) {
                const auto& pair = candidates.pairs[pi];
                const auto left = std::find_if(candidates.left_oligos.begin(), candidates.left_oligos.end(),
                    [&](const Oligo& o) { return same_oligo(o, pair.left); });
                const auto right = std::find_if(candidates.right_oligos.begin(), candidates.right_oligos.end(),
                    [&](const Oligo& o) { return same_oligo(o, pair.right); });
                if (left == candidates.left_oligos.end() || right == candidates.right_oligos.end())
                    throw std::runtime_error("Off-target candidate pair refers to missing primer");
                definitions.push_back({"OUT" + std::to_string(oi) + "_L" +
                    std::to_string(left - candidates.left_oligos.begin()),
                    "OUT" + std::to_string(oi) + "_R" +
                    std::to_string(right - candidates.right_oligos.begin()), false, true, oi, pi});
            }
        }
        const auto products = ac.search_amplicons(ctx.args.ref_file, definitions,
            ctx.args.threshold, ctx.args.dg_thres, ctx.args.off_target_amp_min,
            ctx.args.off_target_amp_max, ctx.args.chunk_size, ctx.args.block_size, ctx.args.nthreads);
        std::vector<std::set<std::size_t>> rejected(all_candidates.size());
        const auto path = ctx.args.output_file + ".offtarget-amplicons.tsv";
        std::ofstream report(path);
        if (!report) throw std::runtime_error("Cannot write " + path);
        report << "contig\tproduct_start_0based\tproduct_end_exclusive\tproduct_bp\toutput_index\tpair_index"
               << "\tforward_hit_label\tforward_hit_start\tforward_hit_dg_cal_per_mol"
               << "\treverse_hit_label\treverse_hit_start\treverse_hit_dg_cal_per_mol\n"
               << std::setprecision(std::numeric_limits<double>::max_digits10);
        for (const auto& product : products) {
            const auto& definition = definitions[product.pair_index];
            rejected[definition.output_index].insert(definition.candidate_pair_index);
            report << product.contig << '\t' << product.start << '\t' << product.end
                   << '\t' << product.end-product.start << '\t' << definition.output_index
                   << '\t' << definition.candidate_pair_index << '\t' << product.forward_hit.label
                   << '\t' << product.forward_hit.ref_index << '\t' << product.forward_hit.dg
                   << '\t' << product.reverse_hit.label << '\t' << product.reverse_hit.ref_index
                   << '\t' << product.reverse_hit.dg << '\n';
        }
        report.close();
        if (!report) throw std::runtime_error("Failed writing " + path);
        all_filtered = all_candidates;
        for (std::size_t oi = 0; oi < all_filtered.size(); ++oi) {
            auto& filtered = all_filtered[oi];
            filtered.pairs.clear();
            for (std::size_t pi = 0; pi < all_candidates[oi].pairs.size(); ++pi)
                if (!rejected[oi].count(pi)) filtered.pairs.push_back(all_candidates[oi].pairs[pi]);
            std::cout << "Off-target output " << oi << ": " << all_candidates[oi].pairs.size()
                      << " -> " << filtered.pairs.size() << " pairs ("
                      << rejected[oi].size() << " rejected by paired products)\n";
            if (filtered.pairs.empty())
                throw std::runtime_error("Paired off-target screening left no valid pairs for output " +
                    std::to_string(oi) + "; details: " + path);
            // Keep all surviving pair scores/metadata; remove only orphan primers.
            auto prune = [&](std::vector<Oligo>& oligos, bool right) {
                oligos.erase(std::remove_if(oligos.begin(), oligos.end(), [&](const Oligo& o) {
                    return std::none_of(filtered.pairs.begin(), filtered.pairs.end(),
                        [&](const PrimerResult& pair) { return same_oligo(o, right ? pair.right : pair.left); });
                }), oligos.end());
            };
            prune(filtered.left_oligos, false); prune(filtered.right_oligos, true);
        }
    }

    // split back to segments using offsets
    for (std::size_t s = 0; s < ctx.segments.size(); ++s) {
        std::size_t start = seg_offsets[s];
        std::size_t end   = (s + 1 < ctx.segments.size())
                          ? seg_offsets[s + 1]
                          : all_filtered.size();

        ctx.segments[s].filtered_primers.clear();
        for (std::size_t j = start; j < end; ++j) {
            PrimerOutput po    = all_filtered[j];
            // Preserve original target metadata when splitting segment results.
            po.pdr_left        = all_candidates[j].pdr_left;
            po.pdr_right       = all_candidates[j].pdr_right;
            po.segment_id      = all_candidates[j].segment_id;
            ctx.segments[s].filtered_primers.push_back(po);
        }
    }
}


void DimerStage::run(PipelineContext& ctx) {
    srand(ctx.args.seed);

    Thal::init(std::string(PRIMER3_PATH) + "/src/primer3_config", 
        ctx.args.mv, 
        ctx.args.dv, 
        ctx.args.dntp, 
        ctx.args.dna_conc, 
        ctx.args.temp);

    using Clock = std::chrono::steady_clock;
    using Ms    = std::chrono::milliseconds;

    // Split filtered_primers into two pools by alternating per segment
    // seg_order[segment_id] = running count so we can alternate
    std::unordered_map<std::size_t, int> seg_order;
    std::vector<PrimerOutput> pool1_primers, pool2_primers;
    // track original index so we can reconstruct later
    std::vector<std::size_t>  pool1_idx,    pool2_idx;

    for (std::size_t i = 0; i < ctx.filtered_primers.size(); i++) {
        const PrimerOutput& po = ctx.filtered_primers[i];
        int& cnt = seg_order[po.segment_id];
        if (cnt % 2 == 0) {
            pool1_primers.push_back(po);
            pool1_idx.push_back(i);
        } else {
            pool2_primers.push_back(po);
            pool2_idx.push_back(i);
        }
        cnt++;
    }


    auto solve_pool = [&](const std::vector<PrimerOutput>& primers,
                      int pool_num) -> std::vector<index_t>
    {
        if (primers.empty()) return {};
        KPartiteGraph g(primers, ctx.args.dimer_temp_min.value_or(ctx.args.temp),
            ctx.args.dimer_temp_max.value_or(ctx.args.temp),
            ctx.args.output_file + ".pool" + std::to_string(pool_num) + ".kpartite");

        const std::size_t selected_edges = std::size_t(g.get_K()) * (g.get_K() - 1) / 2;
        const std::size_t allowed_bad = static_cast<std::size_t>(std::floor(
            static_cast<long double>(selected_edges) * ctx.args.dimer_bad_percent / 100.0L));
        std::cout << "Pool " << pool_num << " rank objective: allowed_bad=" << allowed_bad
                  << " / " << selected_edges << " (" << ctx.args.dimer_bad_percent
                  << "%), maximizing rank " << allowed_bad + 1 << " dG\n";
        const double init_temp = 1000.0;
        const std::size_t sa_iter = ctx.args.iter * 5;
        const double cooling = std::exp(std::log(1e-6 / init_temp) / double(sa_iter));
        struct Comparison {
            std::string method;
            std::vector<index_t> solution;
            double sum, worst, bound;
            long ms;
        };
        std::vector<Comparison> results;
        std::istringstream methods(ctx.args.dimer_methods);
        std::string method;
        std::string selected_method = ctx.args.dimer_select;
        if (selected_method.empty()) {
            std::istringstream candidates(ctx.args.dimer_methods);
            std::string candidate;
            selected_method = ctx.args.dimer_methods.substr(0, ctx.args.dimer_methods.find(','));
            while (std::getline(candidates, candidate, ','))
                if (candidate == "top-k-bound") selected_method = candidate;
        }
        while (std::getline(methods, method, ',')) {
            // Reset per method so method ordering does not change search results.
            srand(ctx.args.seed);
            weight_t objective;
            const auto begin = Clock::now();
            auto solution = method == "sum"
                ? g.solve_sum_sa(4000, sa_iter, init_temp, cooling, objective)
                : g.solve_rank_bottleneck_sa(4000, sa_iter, init_temp, cooling,
                    method == "min-max" ? 0 : allowed_bad, objective);
            const auto ms = std::chrono::duration_cast<Ms>(Clock::now() - begin).count();
            results.push_back({method, solution, g.dimer_sum_score(solution),
                g.rank_bottleneck_score(solution, 0),
                g.rank_bottleneck_score(solution, allowed_bad), ms});
            // Persist each alternative, rather than discarding non-selected panels.
            const std::string path = ctx.args.output_file + ".pool" +
                std::to_string(pool_num) + "." + method + ".primers.tsv";
            std::ofstream file(path);
            if (!file) throw std::runtime_error("Cannot write " + path);
            file << "pool\ttarget_index\tside\tcandidate_index\tsequence\n";
            for (std::size_t i = 0; i < primers.size(); ++i) {
                const auto left = solution[2*i], right = solution[2*i+1];
                file << pool_num << '\t' << i << "\tleft\t" << left << '\t'
                     << primers[i].left_oligos[left].seq << '\n';
                file << pool_num << '\t' << i << "\tright\t" << right << '\t'
                     << primers[i].right_oligos[right].seq << '\n';
            }
            file.close();
            if (!file) throw std::runtime_error("Failed writing " + path);
        }
        const std::string report_path = ctx.args.output_file + ".pool" +
            std::to_string(pool_num) + ".dimer-methods.tsv";
        std::ofstream report(report_path);
        if (!report) throw std::runtime_error("Cannot write " + report_path);
        auto print_table = [&](std::ostream& out) {
            out << "method\tsum_dg_cal_per_mol\tworst_dg_cal_per_mol"
                << "\ttop_k_bound_cal_per_mol\tallowed_bad\tedges\ttime_ms\tselected\n";
            out << std::setprecision(std::numeric_limits<double>::max_digits10);
            for (const auto& r : results)
                out << r.method << '\t' << r.sum << '\t' << r.worst << '\t'
                    << r.bound << '\t' << allowed_bad << '\t' << selected_edges
                    << '\t' << r.ms << '\t' << (r.method == selected_method) << '\n';
        };
        print_table(report);
        report.close();
        if (!report) throw std::runtime_error("Failed writing " + report_path);
        std::cout << "\nPool " << pool_num << " method comparison (higher dG is better)\n";
        std::cout << "dG: cal/mol; time: ms\n";
        auto format_score = [](double value) {
            std::ostringstream text;
            text << std::fixed << std::setprecision(2) << value;
            return text.str();
        };
        std::vector<std::vector<std::string>> rows = {
            {"Method", "Sum dG", "Worst dG", "Top-k bound", "Allowed bad", "Edges", "Time(ms)", "Selected"}
        };
        for (const auto& r : results)
            rows.push_back({r.method, format_score(r.sum), format_score(r.worst),
                format_score(r.bound), std::to_string(allowed_bad),
                std::to_string(selected_edges), std::to_string(r.ms),
                r.method == selected_method ? "yes" : ""});
        std::vector<std::size_t> widths(rows.front().size(), 0);
        for (const auto& row : rows)
            for (std::size_t col = 0; col < row.size(); ++col)
                widths[col] = std::max(widths[col], row[col].size());
        auto print_console_row = [&](const std::vector<std::string>& row) {
            std::ostringstream line;
            for (std::size_t col = 0; col < row.size(); ++col) {
                if (col) line << "  ";
                line << (col == 0 || col == row.size() - 1 ? std::left : std::right)
                     << std::setw(static_cast<int>(widths[col])) << row[col];
            }
            std::cout << line.str() << '\n';
        };
        print_console_row(rows.front());
        std::vector<std::string> divider;
        for (auto width : widths) divider.push_back(std::string(width, '-'));
        print_console_row(divider);
        for (std::size_t row = 1; row < rows.size(); ++row)
            print_console_row(rows[row]);
        for (const auto& r : results)
            if (r.method == selected_method) {
                std::cout << " Selected: " << r.method << '\n';
                return r.solution;
            }
        throw std::runtime_error("No selected dimer method result");
    };

    auto sol1 = solve_pool(pool1_primers, 1);
    auto sol2 = solve_pool(pool2_primers, 2);

    // Reconstruct PrimerResult for pool 1
    for (std::size_t amp = 0; amp < pool1_primers.size(); amp++) {
        index_t left_n  = sol1[amp * 2];
        index_t right_n = sol1[amp * 2 + 1];

        const PrimerOutput& po = pool1_primers[amp];
        PrimerResult result;
        result.left         = po.left_oligos[left_n];
        result.right        = po.right_oligos[right_n];
        const auto selected_pair = std::find_if(po.pairs.begin(), po.pairs.end(),
            [&](const PrimerResult& pair) {
                return same_oligo(pair.left, result.left) && same_oligo(pair.right, result.right);
            });
        if (po.pairs_required && selected_pair == po.pairs.end())
            throw std::runtime_error("Dimer result does not belong to retained pair candidates");
        if (selected_pair != po.pairs.end()) result.msa = selected_pair->msa;
        result.product_size = result.right.start - result.left.start + result.right.length;
        result.pdr_left     = po.pdr_left;
        result.pdr_right    = po.pdr_right;
        result.segment_id   = po.segment_id;
        result.pool         = 1;

        ctx.solution_primers.push_back(result);
    }

    // Reconstruct PrimerResult for pool 2
    for (std::size_t amp = 0; amp < pool2_primers.size(); amp++) {
        index_t left_n  = sol2[amp * 2];
        index_t right_n = sol2[amp * 2 + 1];

        const PrimerOutput& po = pool2_primers[amp];
        PrimerResult result;
        result.left         = po.left_oligos[left_n];
        result.right        = po.right_oligos[right_n];
        const auto selected_pair = std::find_if(po.pairs.begin(), po.pairs.end(),
            [&](const PrimerResult& pair) {
                return same_oligo(pair.left, result.left) && same_oligo(pair.right, result.right);
            });
        if (po.pairs_required && selected_pair == po.pairs.end())
            throw std::runtime_error("Dimer result does not belong to retained pair candidates");
        if (selected_pair != po.pairs.end()) result.msa = selected_pair->msa;
        result.product_size = result.right.start - result.left.start + result.right.length;
        result.pdr_left     = po.pdr_left;
        result.pdr_right    = po.pdr_right;
        result.segment_id   = po.segment_id;
        result.pool         = 2;

        ctx.solution_primers.push_back(result);
    }

    // Sort solution_primers back to original amplicon order
    std::vector<std::size_t> all_idx;
    all_idx.insert(all_idx.end(), pool1_idx.begin(), pool1_idx.end());
    all_idx.insert(all_idx.end(), pool2_idx.begin(), pool2_idx.end());

    std::vector<PrimerResult> sorted(ctx.solution_primers.size());
    for (std::size_t i = 0; i < all_idx.size(); i++)
        sorted[all_idx[i]] = ctx.solution_primers[i];
    ctx.solution_primers = std::move(sorted);
}

void Pipeline::run(PipelineContext& ctx) {
    constexpr std::size_t W = 80;

    auto shorten_name = [](const std::string& name,
                           std::size_t max_len) -> std::string {
        if (name.size() <= max_len) {
            return name;
        }

        constexpr std::string_view separator = "...";

        if (max_len <= separator.size()) {
            return name.substr(0, max_len);
        }

        const std::size_t available = max_len - separator.size();
        const std::size_t prefix_len = (available + 1) / 2;
        const std::size_t suffix_len = available / 2;

        return name.substr(0, prefix_len) +
               std::string(separator) +
               name.substr(name.size() - suffix_len);
    };

    auto print_banner = [W](const std::string& text) {
        const std::size_t width = std::max(W, text.size());
        const std::size_t left_pad = (width - text.size()) / 2;
        const std::size_t right_pad =
            width - text.size() - left_pad;

        std::cout << std::string(left_pad, '=')
                  << text
                  << std::string(right_pad, '=')
                  << '\n';
    };

    std::cout << '\n'
              << std::string(W, '=') << '\n'
              << std::string((W - 26) / 2, ' ')
              << "PIPELINE EXECUTION STARTED\n"
              << std::string(W, '=') << '\n';

    for (std::size_t i = 0; i < stages_.size(); ++i) {
        const std::string stage_number =
            std::to_string(i + 1) + "/" +
            std::to_string(stages_.size());

        const std::string short_name = stages_[i]->short_name();

        if (short_name == "dim") {
            // Flatten all filtered primers from all segments.
            ctx.filtered_primers.clear();

            for (const auto& seg : ctx.segments) {
                ctx.filtered_primers.insert(
                    ctx.filtered_primers.end(),
                    seg.filtered_primers.begin(),
                    seg.filtered_primers.end()
                );
            }

            const std::string running =
                "[ " + stage_number +
                " Running: " + stages_[i]->name() +
                " (global) ]";

            const std::string completed =
                "[ Completed: " + stages_[i]->name() + " ]";

            std::cout << "\n\n";
            print_banner(running);

            stages_[i]->run(ctx);

            print_banner(completed);
            std::cout << '\n';

        } else if (short_name == "offt") {
            // Global stage operating directly on ctx.segments.
            const std::string running =
                "[ " + stage_number +
                " Running: " + stages_[i]->name() +
                " (global) ]";

            const std::string completed =
                "[ Completed: " + stages_[i]->name() + " ]";

            std::cout << "\n\n";
            print_banner(running);

            stages_[i]->run(ctx);

            print_banner(completed);
            std::cout << '\n';

        } else {
            for (std::size_t s = 0; s < ctx.segments.size(); ++s) {
                ctx.load_segment(s);

                const std::string segment_counter =
                    " [" + std::to_string(s + 1) +
                    "/" + std::to_string(ctx.segments.size()) + "]";

                const std::string running_prefix =
                    "[ " + stage_number +
                    " Running: " + stages_[i]->name() +
                    " | ";

                const std::string running_suffix =
                    segment_counter + " ]";

                // Restrict the segment name so the running banner fits in W.
                const std::size_t fixed_running_length =
                    running_prefix.size() + running_suffix.size();

                const std::size_t max_segment_name_length =
                    fixed_running_length < W
                        ? W - fixed_running_length
                        : 15;

                const std::string segment_name =
                    shorten_name(
                        ctx.segments[s].name,
                        max_segment_name_length
                    );

                const std::string segment_tag =
                    segment_name + segment_counter;

                const std::string running =
                    running_prefix + segment_tag + " ]";

                const std::string completed_prefix =
                    "[ Completed: " + stages_[i]->name() +
                    " | ";

                const std::string completed =
                    completed_prefix + segment_tag + " ]";

                std::cout << "\n\n";
                print_banner(running);

                stages_[i]->run(ctx);
                ctx.save_segment();

                print_banner(completed);
                std::cout << '\n';
            }
        }
    }
}
