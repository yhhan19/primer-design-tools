#include "pipeline.hpp"
#include "pdr_layout.hpp"
#include "run_log.hpp"

namespace {
using Interval = std::pair<int, int>;
Interval span(const PrimerResult& p) { return {p.left.start, p.right.start + p.right.length}; }
std::vector<Interval> merge(std::vector<Interval> intervals) {
    std::sort(intervals.begin(), intervals.end());
    std::vector<Interval> result;
    for (auto interval : intervals) {
        if (!result.empty() && interval.first <= result.back().second)
            result.back().second = std::max(result.back().second, interval.second);
        else result.push_back(interval);
    }
    return result;
}
void add_oligo(std::vector<Oligo>& oligos, const Oligo& o) {
    if (std::none_of(oligos.begin(), oligos.end(), [&](const Oligo& v) { return same_oligo(v, o); }))
        oligos.push_back(o);
}
void add_pair(PrimerOutput& output, const PrimerResult& p) {
    if (std::none_of(output.pairs.begin(), output.pairs.end(), [&](const PrimerResult& v) {
        return same_oligo(v.left, p.left) && same_oligo(v.right, p.right);
    })) {
        output.pairs.push_back(p);
        add_oligo(output.left_oligos, p.left);
        add_oligo(output.right_oligos, p.right);
    }
}
std::vector<double> objective(const std::vector<PrimerResult>& panel, const Args& args, int pool) {
    std::vector<const PrimerResult*> selected;
    for (const auto& p : panel) if (p.pool == std::size_t(pool)) selected.push_back(&p);
    if (selected.empty()) return {0};
    std::vector<std::string> primers;
    for (const auto* p : selected) {
        primers.push_back(p->left.seq); primers.push_back(revcomp(p->right.seq));
    }
    std::vector<double> edges;
    std::vector<bool> affected(primers.size(), false);
    const auto cutoff = args.dimer_dg_cutoff.value_or(Args::default_dimer_dg_cutoff);
    double severity = 0, sum = 0, bad = 0;
    for (std::size_t p = 0; p < primers.size(); ++p)
        for (std::size_t q = p + 1; q < primers.size(); ++q) {
            const double dg = Thal::compute_dimer_worst_dg(primers[p], primers[q],
                args.dimer_temp_min.value_or(args.temp), args.dimer_temp_max.value_or(args.temp));
            edges.push_back(dg); sum += dg;
            if (dg < cutoff) { affected[p] = affected[q] = true; ++bad; severity += cutoff - dg; }
        }
    std::string method = args.dimer_select;
    if (method.empty()) method = args.dimer_methods.find("top-k-bound") != std::string::npos
        ? "top-k-bound" : args.dimer_methods.substr(0, args.dimer_methods.find(','));
    if (method == "min-affected-primers")
        return {-double(std::count(affected.begin(), affected.end(), true))};
    if (method == "min-affected-amplicons") {
        double count = 0;
        for (std::size_t i = 0; i < affected.size(); i += 2) count += affected[i] || affected[i+1];
        return {-count, -severity};
    }
    if (method == "min-below-cutoff") return {-bad};
    if (method == "sum") return {sum};
    const std::size_t allowed = method == "top-k-bound"
        ? std::size_t(std::floor(static_cast<long double>(edges.size()) * args.dimer_bad_percent / 100)) : 0;
    std::nth_element(edges.begin(), edges.begin() + allowed, edges.end());
    return {edges[allowed]};
}
}

bool expansion_layout_compatible(const std::vector<PrimerResult>& original,
                                 const std::vector<PrimerResult>& trial) {
    if (original.size() != trial.size()) return false;
    std::map<std::size_t, std::vector<Interval>> before, after;
    std::map<std::pair<std::size_t, std::size_t>, std::vector<Interval>> pool_spans;
    for (std::size_t i = 0; i < original.size(); ++i) {
        if (original[i].segment_id != trial[i].segment_id || original[i].pool != trial[i].pool)
            return false;
        before[original[i].segment_id].push_back(span(original[i]));
        after[trial[i].segment_id].push_back(span(trial[i]));
        pool_spans[{trial[i].segment_id, trial[i].pool}].push_back(span(trial[i]));
        if (trial[i].left.start < 0 || trial[i].left.length <= 0 ||
            trial[i].right.length <= 0 || span(trial[i]).second <= trial[i].left.start)
            return false;
    }
    for (auto& entry : pool_spans) {
        auto& intervals = entry.second;
        std::sort(intervals.begin(), intervals.end());
        for (std::size_t i = 1; i < intervals.size(); ++i)
            if (intervals[i-1].second > intervals[i].first) return false;
    }
    for (const auto& segment : before) {
        const auto coverage = merge(after[segment.first]);
        for (const auto& interval : merge(segment.second))
            if (std::none_of(coverage.begin(), coverage.end(), [&](Interval v) {
                return v.first <= interval.first && v.second >= interval.second;
            })) return false;
    }
    return true;
}

void PDRExpansionStage::run(PipelineContext& ctx) {
    if (ctx.args.pdr_search != "midpoint") throw std::runtime_error("Expansion requires midpoint mode");
    if (ctx.solution_primers.size() != ctx.filtered_primers.size())
        throw std::runtime_error("Expansion requires a complete current panel");
    const auto prefix = ctx.args.output_file;
    auto original = ctx.solution_primers;
    annotate_solution_dimers(ctx.segments, original,
        ctx.args.dimer_dg_cutoff.value_or(Args::default_dimer_dg_cutoff),
        ctx.args.dimer_temp_min.value_or(ctx.args.temp),
        ctx.args.dimer_temp_max.value_or(ctx.args.temp), prefix + ".before-expansion");
    save_solution(ctx.segments, original, prefix + ".before-expansion.bed");
    save_coverage_report(ctx.segments, original, prefix + ".before-expansion");
    std::ofstream report(prefix + ".expansion.tsv");
    if (!report) throw std::runtime_error("Cannot write expansion audit");
    report << "target\tsegment\tleft_start\tleft_end\tright_start\tright_end\tnew_pairs_before_screening\n";
    PipelineContext screening;
    screening.args = ctx.args;
    screening.args.output_file = prefix + ".expanded";
    screening.segments = ctx.segments;
    for (auto& s : screening.segments) { s.candidate_primers.clear(); s.filtered_primers.clear(); }
    std::vector<std::vector<std::size_t>> mapping(ctx.segments.size());
    std::vector<std::vector<PDRLayoutWindow>> layouts(ctx.segments.size());
    std::vector<std::size_t> layout_index(original.size());
    for (std::size_t i = 0; i < original.size(); ++i) {
        const auto& po = ctx.filtered_primers[i];
        layout_index[i] = layouts.at(po.segment_id).size();
        layouts[po.segment_id].push_back({po.pdr_left, std::int64_t(po.pdr_left) + std::int64_t(ctx.args.len_PDR),
                                         po.pdr_right, std::int64_t(po.pdr_right) + std::int64_t(ctx.args.len_PDR)});
    }
    for (std::size_t sid = 0; sid < layouts.size(); ++sid) {
        const auto error = expanded_pdr_layout_error(layouts[sid], ctx.segments[sid].tmpl.size(),
                                             ctx.args.len_fu, ctx.args.len_ru);
        if (!error.empty()) throw std::runtime_error("Expansion source layout: " + error);
    }
    for (std::size_t i = 0; i < original.size(); ++i)
        expand_pdr_window(layouts[original[i].segment_id], layout_index[i],
            !original[i].left.dimer_risk.partners.empty(), !original[i].right.dimer_risk.partners.empty(),
            ctx.args.pdr_expansion_bases, ctx.segments[original[i].segment_id].tmpl.size(),
            ctx.args.len_fu, ctx.args.len_ru);
    // Export every window, including unchanged targets and rejected trials.
    std::ofstream window_audit(prefix + ".expanded.pdr-windows.tsv");
    if (!window_audit) throw std::runtime_error("Cannot write trial PDR windows");
    window_audit << "segment\ttarget\tleft_start\tleft_end\tright_start\tright_end\n";
    for (std::size_t sid = 0; sid < layouts.size(); ++sid)
        for (std::size_t n = 0; n < layouts[sid].size(); ++n) {
            const auto& w = layouts[sid][n];
            window_audit << sid << '\t' << n << '\t' << w.left_start << '\t' << w.left_end
                         << '\t' << w.right_start << '\t' << w.right_end << '\n';
        }
    window_audit.close();
    if (!window_audit) throw std::runtime_error("Failed writing trial PDR windows");
    for (std::size_t i = 0; i < original.size(); ++i) {
        const bool left_bad = !original[i].left.dimer_risk.partners.empty();
        const bool right_bad = !original[i].right.dimer_risk.partners.empty();
        if (!left_bad && !right_bad) continue;
        const auto& old = ctx.filtered_primers[i];
        const auto sid = old.segment_id;
        const auto& expanded = layouts[sid][layout_index[i]];
        const auto l = std::make_pair(std::size_t(expanded.left_start), std::size_t(expanded.left_end));
        const auto r = std::make_pair(std::size_t(expanded.right_start), std::size_t(expanded.right_end));
        if (l.first == old.pdr_left && l.second == old.pdr_left + ctx.args.len_PDR &&
            r.first == old.pdr_right && r.second == old.pdr_right + ctx.args.len_PDR) continue;
        PipelineContext generation;
        RunLogBlock generation_block("Expand " + ctx.segments[sid].name + " / amplicon " + std::to_string(layout_index[i]+1));
        generation.args = screening.args;
        generation.args.output_file += ".target" + std::to_string(i);
        generation.segments.resize(sid + 1);
        generation.load_segment(sid);
        generation.tmpl = ctx.segments[sid].tmpl;
        generation.labels = ctx.segments[sid].labels;
        generation.sequences = ctx.segments[sid].sequences;
        generation.pdr_regions = {index_t(l.first), index_t(r.first)};
        generation.pdr_region_lengths = {l.second - l.first, r.second - r.first};
        PrimerSelStage generator;
        generator.run(generation);
        PrimerOutput combined = old;
        const auto before = combined.pairs.size();
        for (const auto& pair : generation.candidate_primers.at(0).pairs) {
            const auto length = pair.right.start + pair.right.length - pair.left.start;
            if (length < int(ctx.args.len_amp_min) || length > int(ctx.args.len_amp)) continue;
            auto replacement = original;
            replacement[i].left = pair.left; replacement[i].right = pair.right;
            if (expansion_layout_compatible(original, replacement)) add_pair(combined, pair);
        }
        report << i << '\t' << sid << '\t' << l.first << '\t' << l.second << '\t'
               << r.first << '\t' << r.second << '\t' << combined.pairs.size() - before << '\n';
        if (combined.pairs.size() == before) {
            generation_block.complete("no additional compatible pairs");
            continue;
        }
        generation_block.complete(std::to_string(combined.pairs.size()-before) + " new compatible pairs");
        screening.segments[sid].candidate_primers.push_back(combined);
        mapping[sid].push_back(i);
    }
    bool generated = false;
    for (const auto& map : mapping) generated |= !map.empty();
    if (!generated) { report << "# status\tno_new_compatible_pairs\n"; return; }
    OffTargetStage filter;
    filter.run(screening);
    PipelineContext trial = ctx;
    trial.args.output_file = prefix + ".expanded";
    trial.args.dimer_start_prefix = prefix; // Old candidate indices are kept stable below.
    trial.solution_primers.clear();
    std::size_t added = 0;
    for (std::size_t sid = 0; sid < mapping.size(); ++sid)
        for (std::size_t j = 0; j < mapping[sid].size(); ++j) {
            auto& output = trial.filtered_primers[mapping[sid][j]];
            const auto before = output.pairs.size();
            for (const auto& p : screening.segments[sid].filtered_primers.at(j).pairs) add_pair(output, p);
            added += output.pairs.size() - before;
        }
    if (!added) { report << "# status\tno_new_pairs_after_screening\n"; return; }
    for (std::size_t i = 0; i < trial.filtered_primers.size(); ++i) {
        auto& output = trial.filtered_primers[i];
        const auto& w = layouts[original[i].segment_id][layout_index[i]];
        output.pdr_left = w.left_start;
        output.pdr_right = w.right_start;
        const auto before_count = output.pairs.size();
        // Apply constraints to all old and new pair choices. Keep oligo arrays
        // stable so incumbent graph-node indices remain valid; only allowed
        // pair choices are restricted. In particular, terminal candidates
        // cannot lose bases covered by the original panel.
        output.pairs.erase(std::remove_if(output.pairs.begin(), output.pairs.end(),
            [&](const PrimerResult& p) {
                const auto product = span(p);
                if (p.left.start < w.left_start || p.left.start + p.left.length > w.left_end ||
                    p.right.start < w.right_start || p.right.start + p.right.length > w.right_end ||
                    product.second - product.first < int(ctx.args.len_amp_min) ||
                    product.second - product.first > int(ctx.args.len_amp)) return true;
                auto replacement = original;
                replacement[i].left = p.left; replacement[i].right = p.right;
                return !expansion_layout_compatible(original, replacement);
            }), output.pairs.end());
        report << "# constraint_filtered_pairs\t" << i << '\t' << before_count - output.pairs.size() << '\n';
        if (output.pairs.empty()) throw std::runtime_error("Expansion removed the incumbent candidate pair");
        output.pairs_required = true;
    }
    DimerStage selector;
    selector.run(trial);
    annotate_solution_dimers(trial.segments, trial.solution_primers,
        ctx.args.dimer_dg_cutoff.value_or(Args::default_dimer_dg_cutoff),
        ctx.args.dimer_temp_min.value_or(ctx.args.temp),
        ctx.args.dimer_temp_max.value_or(ctx.args.temp), trial.args.output_file);
    save_solution(trial.segments, trial.solution_primers, trial.args.output_file + ".bed");
    save_coverage_report(trial.segments, trial.solution_primers, trial.args.output_file);
    bool improved = false, nonworse = true;
    for (int pool : {1, 2}) {
        const auto before = objective(original, ctx.args, pool);
        const auto after = objective(trial.solution_primers, ctx.args, pool);
        improved |= after > before; nonworse &= after >= before;
    }
    bool valid = expansion_layout_compatible(original, trial.solution_primers);
    std::string layout_reason = valid ? "" : "lost covered template positions, same-pool overlap or target identity";
    for (const auto& p : trial.solution_primers) {
        const auto interval = span(p);
        const auto length = interval.second - interval.first;
        if (!(length >= int(ctx.args.len_amp_min) && length <= int(ctx.args.len_amp) &&
              interval.first >= 0 && interval.second <= int(ctx.segments[p.segment_id].tmpl.size()))) {
            valid = false; layout_reason = "selected product length or target boundary";
        }
    }
    for (std::size_t sid = 0; sid < layouts.size(); ++sid) {
        const auto error = expanded_pdr_layout_error(layouts[sid], ctx.segments[sid].tmpl.size(),
                                            ctx.args.len_fu, ctx.args.len_ru);
        if (!error.empty()) { valid = false; layout_reason = error; }
    }
    for (std::size_t i = 0; i < trial.solution_primers.size(); ++i) {
        const auto& p = trial.solution_primers[i];
        const auto& w = layouts[p.segment_id][layout_index[i]];
        if (!(p.left.start >= w.left_start && p.left.start + p.left.length <= w.left_end &&
              p.right.start >= w.right_start && p.right.start + p.right.length <= w.right_end)) {
            valid = false; layout_reason = "selected primer outside validated expanded window";
        }
    }
    if (!valid || !nonworse || !improved) {
        report << "# status\t" << (!valid ? "rejected_combined_layout" : "retained_original_objective") << '\n';
        if (!valid) report << "# layout_reason\t" << layout_reason << '\n';
        return;
    }
    // Preserve old graph/method exports and promote the new graph family.
    const auto directory = fs::path(prefix).parent_path().empty() ? fs::path(".") : fs::path(prefix).parent_path();
    const auto old_name = fs::path(prefix).filename().string();
    std::vector<fs::path> files;
    for (const auto& entry : fs::directory_iterator(directory))
        if (entry.is_regular_file() && entry.path().filename().string().find(old_name + ".expanded.pool") == 0)
            files.push_back(entry.path());
    for (const auto& source : files) {
        const auto tail = source.filename().string().substr((old_name + ".expanded").size());
        if (!ctx.args.dimer_landscape_points && tail.find(".landscape.") != std::string::npos) continue;
        const auto target = fs::path(prefix + tail);
        if (fs::exists(target))
            fs::copy_file(target, prefix + ".before-expansion" + tail, fs::copy_options::overwrite_existing);
        fs::copy_file(source, target, fs::copy_options::overwrite_existing);
    }
    ctx.filtered_primers = std::move(trial.filtered_primers);
    if (ctx.args.dimer_landscape_points && fs::exists(prefix + ".expanded.landscape.svg")) {
        if (fs::exists(prefix + ".landscape.svg"))
            fs::copy_file(prefix + ".landscape.svg", prefix + ".before-expansion.landscape.svg", fs::copy_options::overwrite_existing);
        fs::copy_file(prefix + ".expanded.landscape.svg", prefix + ".landscape.svg", fs::copy_options::overwrite_existing);
    }
    ctx.solution_primers = std::move(trial.solution_primers);
    std::ofstream windows(prefix + ".pdr-windows.tsv");
    windows << "segment\ttarget\tleft_start\tleft_end\tright_start\tright_end\n";
    for (std::size_t sid = 0; sid < layouts.size(); ++sid) {
        ctx.segments[sid].pdr_regions.clear();
        ctx.segments[sid].pdr_region_lengths.clear();
        for (std::size_t i = 0; i < layouts[sid].size(); ++i) {
            const auto& w = layouts[sid][i];
            ctx.segments[sid].pdr_regions.insert(ctx.segments[sid].pdr_regions.end(),
                {index_t(w.left_start), index_t(w.right_start)});
            ctx.segments[sid].pdr_region_lengths.insert(ctx.segments[sid].pdr_region_lengths.end(),
                {std::size_t(w.left_end - w.left_start), std::size_t(w.right_end - w.right_start)});
            windows << sid << '\t' << i << '\t' << w.left_start << '\t' << w.left_end
                    << '\t' << w.right_start << '\t' << w.right_end << '\n';
        }
    }
    windows.close();
    if (!windows) throw std::runtime_error("Failed writing expanded PDR windows");
    for (auto& s : ctx.segments) s.filtered_primers.clear();
    std::map<std::size_t, std::size_t> ordinal;
    for (const auto& po : ctx.filtered_primers) {
        const auto n = ordinal[po.segment_id]++;
        ctx.segments[po.segment_id].filtered_primers.push_back(po);
        const auto path = prefix + ".segment" + std::to_string(po.segment_id) + ".pdr" + std::to_string(n) + ".pair-msa.tsv";
        if (fs::exists(path)) fs::copy_file(path, path + ".before-expansion", fs::copy_options::overwrite_existing);
        std::ofstream msa(path);
        if (!msa) throw std::runtime_error("Cannot write expanded MSA scores");
        msa << "rank\tleft_start\tright_start\tleft_reference_sequence\tright_reference_sequence"
            << "\tjoint_coverage\tjoint_perfect_match\tjoint_3prime_match\tweighted_mismatch\tunknown_fraction\tmsa_sequences\tkept\n";
        std::size_t rank = 0;
        for (const auto& pair : po.pairs) {
            const auto& m = pair.msa;
            msa << std::setprecision(17) << ++rank << '\t' << pair.left.start << '\t' << pair.right.start
                << '\t' << pair.left.seq << '\t' << pair.right.seq << '\t' << m.joint_coverage
                << '\t' << m.joint_perfect_match << '\t' << m.joint_3prime_match << '\t'
                << m.weighted_mismatch << '\t' << m.unknown_fraction << '\t' << m.sequence_count << "\t1\n";
        }
    }
    for (std::size_t sid = 0; sid < ctx.segments.size(); ++sid)
        for (std::size_t n = 0; n < ctx.segments[sid].filtered_primers.size(); ++n) {
            const auto path = prefix + ".segment" + std::to_string(sid) +
                ".pdr" + std::to_string(n) + ".pairs.tsv";
            if (fs::exists(path)) fs::copy_file(path, path + ".before-expansion", fs::copy_options::overwrite_existing);
            std::ofstream pairs(path);
            pairs << "pair_index\tleft_start\tright_start\tleft_sequence\tright_reference_sequence\tproduct_bp\n";
            std::size_t index = 0;
            for (const auto& p : ctx.segments[sid].filtered_primers[n].pairs)
                pairs << index++ << '\t' << p.left.start << '\t' << p.right.start << '\t'
                      << p.left.seq << '\t' << p.right.seq << '\t'
                      << p.right.start + p.right.length - p.left.start << '\n';
            pairs.close();
            if (!pairs) throw std::runtime_error("Failed writing expanded pair list");
        }
    report << "# status\taccepted\n# added_filtered_pairs\t" << added << '\n';
}
