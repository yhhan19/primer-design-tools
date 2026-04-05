#include "pipeline.hpp"

void PDRStage::run(PipelineContext& ctx) {
    srand(ctx.args.seed);

    RiskOptimizer ro(ctx.labels,
                     ctx.sequences, 
                     ctx.args.len_PDR, 
                     ctx.args.len_amp, 
                     ctx.args.len_amp);

    auto PDR = ro.search(0, ctx.args.u_max, ctx.args.u_min);

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

    // --- 1. Global settings ---
    p3_global_settings *pa = p3_create_global_settings();
    p3_set_gs_primer_opt_size(pa, 20);
    p3_set_gs_primer_min_size(pa, 18);
    p3_set_gs_primer_max_size(pa, 25);
    p3_set_gs_primer_opt_tm(pa, 60.0);
    p3_set_gs_primer_min_tm(pa, 57.0);
    p3_set_gs_primer_max_tm(pa, 63.0);
    p3_set_gs_primer_min_gc(pa, 40.0);
    p3_set_gs_primer_max_gc(pa, 60.0);

    pa->p_args.salt_conc     = 50.0;
    pa->p_args.divalent_conc = 1.5;
    pa->p_args.dntp_conc     = 0.6;
    pa->p_args.dna_conc      = 50.0;

    pa->pr_min[0] = ctx.args.len_amp_min;
    pa->pr_max[0] = ctx.args.len_amp + ctx.args.len_PDR;
    pa->num_intervals = 1;

    pa->pick_left_primer    = 1;
    pa->pick_right_primer   = 1;
    pa->pick_internal_oligo = 0;
    pa->num_return          = 10;

    std::cout << "Product size: [" << pa->pr_min[0] 
              << ", " << pa->pr_max[0] << "] bp"
              << "  |  Template length: " << ctx.tmpl.size() << " bp\n";

    // --- 2. Sequence arguments ---
    seq_args *sa = create_seq_arg();
    sa->sequence      = strdup(ctx.tmpl.c_str());
    sa->sequence_name = strdup("my_template");

    sa->ok_regions.any_left  = 0;
    sa->ok_regions.any_right = 0;
    sa->ok_regions.any_pair  = 0;
    sa->ok_regions.count     = 1;
    sa->ok_regions.left_pairs[0][1]  = ctx.args.len_PDR;
    sa->ok_regions.right_pairs[0][1] = ctx.args.len_PDR;

    // --- 3. Run primer design ---
    const int TW = 10 + 20 + 10 + 10 + 10; // = 60
    std::cout << std::string(TW, '-') << "\n"
              << std::right
              << std::setw(10) << "Index"
              << std::setw(20) << "PDRs"
              << std::setw(10) << "Left"
              << std::setw(10) << "Right"
              << std::setw(10) << "Pairs"
              << "\n"
              << std::string(TW, '-') << "\n";

    for (std::size_t i = 0; i < ctx.pdr_regions.size(); i += 2) {
        int region_idx = i / 2;
        sa->ok_regions.left_pairs[0][0]  = ctx.pdr_regions[i];
        sa->ok_regions.right_pairs[0][0] = ctx.pdr_regions[i + 1];

        p3retval *retval = choose_primers(pa, sa);

        if (retval == nullptr) {
            std::cerr << "ERROR: choose_primers() returned NULL for region "
                      << region_idx << "\n";
            return;
        }
        if (retval->glob_err.data)
            std::cerr << "Global error: "   
                      << retval->glob_err.data << "\n";
        if (retval->per_sequence_err.data)
            std::cerr << "Sequence error: " 
                      << retval->per_sequence_err.data << "\n";

        int n_pairs = retval->best_pairs.num_pairs;
        std::string range = "[" + std::to_string(ctx.pdr_regions[i]) +
                            ", " + std::to_string(ctx.pdr_regions[i+1]) + "]";
        std::cout << std::right
                  << std::setw(10) << region_idx
                  << std::setw(20) << range
                  << std::setw(10) << retval->fwd.num_elem
                  << std::setw(10) << retval->rev.num_elem
                  << std::setw(10) << n_pairs
                  << "\n";

        ctx.candidate_primers.push_back(extract_all(retval, sa));
        destroy_p3retval(retval);
    }

    std::cout << std::string(TW, '-') << "\n"
              << "Total PDRs processed: " 
              << ctx.pdr_regions.size() / 2 << "\n";

    destroy_seq_args(sa);
    p3_destroy_global_settings(pa);
}

void DimerStage::run(PipelineContext& ctx) {
    srand(ctx.args.seed);

    Thal::init(std::string(PRIMER3_PATH) + "/src/primer3_config", 
        ctx.args.mv, 
        ctx.args.dv, 
        ctx.args.dntp, 
        ctx.args.dna_conc, 
        ctx.args.temp);

    KPartiteGraph g(ctx.filtered_primers);
    auto solution = g.solve_fast(ctx.args.iter);
    std::cout << "loss: " << g.cost(solution) << std::endl;

    ctx.dimer_solution = solution;
}

void OffTargetStage::run(PipelineContext& ctx) {
    if (ctx.args.ref_file == "") {
        std::cout << "No ref file provided, skip off target search\n";
        ctx.filtered_primers = ctx.candidate_primers;
        return ;
    }

    Thal::init(std::string(PRIMER3_PATH) + "/src/primer3_config", 
            ctx.args.mv, 
            ctx.args.dv, 
            ctx.args.dntp, 
            ctx.args.dna_conc, 
            ctx.args.temp);
    
    std::vector<std::string> labels, sequences;
    convert(ctx.candidate_primers, labels, sequences);

    Automaton *ac = new Automaton(labels, sequences, ctx.args.kmer_len);
    auto results = ac->search(ctx.args.ref_file, 
            ctx.args.kmer_len - 1, 
            ctx.args.threshold,
            ctx.args.dg_thres,
            ctx.args.chunk_size, 
            ctx.args.block_size, 
            ctx.args.nthreads);
    delete ac;

    write_results(ctx.args.output_file + "." + short_name(), results, labels);

    ctx.filtered_primers = filterByDG(results, -18000, ctx.candidate_primers);
    // for (auto &out : ctx.filtered_primers) display_primer_output(out);
}

void Pipeline::run(PipelineContext& ctx) {
    const int W = 80;
    std::cout << "\n"
              << std::string(W, '=') << "\n"
              << std::string((W - 26) / 2, ' ') 
              << "PIPELINE EXECUTION STARTED\n"
              << std::string(W, '=') << "\n";

    for (size_t i = 0; i < stages_.size(); ++i) {
        std::string running   = "[ " + std::to_string(i+1) + 
                                "/" + std::to_string(stages_.size()) + 
                                " Running: " + stages_[i]->name() + " ]";
        std::string completed = "[ Completed: " + stages_[i]->name() + " ]";

        int rpad = (W - running.size())   / 2;
        int cpad = (W - completed.size()) / 2;

        std::cout << "\n\n"
                  << std::string(rpad, '=') << running 
                  << std::string(W - rpad - running.size(),   '=') << "\n";

        stages_[i]->run(ctx);

        std::cout << std::string(cpad, '=') << completed 
                  << std::string(W - cpad - completed.size(), '=') << "\n\n";
    }
}
