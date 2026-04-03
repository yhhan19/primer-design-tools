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

Oligo PrimerSelStage::extract_oligo(const primer_rec& p,
                                    const char* sequence,
                                    bool is_right) {
    Oligo o;
    o.length = p.length;
    o.tm     = p.temp;
    o.gc     = p.gc_content;
    o.start  = is_right ? p.start - p.length + 1 : p.start;
    o.seq    = std::string(sequence + o.start, p.length);
    return o;
}

PrimerOutput PrimerSelStage::extract_all(const p3retval* retval,
                                         const seq_args*  sa) {
    PrimerOutput out;

    out.pairs.reserve(retval->best_pairs.num_pairs);
    for (int i = 0; i < retval->best_pairs.num_pairs; i++) {
        const primer_pair& pp = retval->best_pairs.pairs[i];
        PrimerResult r;
        r.left         = extract_oligo(*pp.left,  sa->sequence, false);
        r.right        = extract_oligo(*pp.right, sa->sequence, true);
        r.product_size = pp.product_size;
        out.pairs.push_back(r);
    }

    out.left_oligos.reserve(retval->fwd.num_elem);
    for (int i = 0; i < retval->fwd.num_elem; i++)
        out.left_oligos.push_back(
            extract_oligo(retval->fwd.oligo[i], sa->sequence, false));

    out.right_oligos.reserve(retval->rev.num_elem);
    for (int i = 0; i < retval->rev.num_elem; i++)
        out.right_oligos.push_back(
            extract_oligo(retval->rev.oligo[i], sa->sequence, true));

    return out;
}

void PrimerSelStage::run(PipelineContext& ctx) {
    // --- 1. Global settings (primer constraints) ---
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

    // Product size range: 100–300 bp
    pa->pr_min[0] = ctx.args.len_amp_min;
    pa->pr_max[0] = ctx.args.len_amp + ctx.args.len_PDR;
    pa->num_intervals = 1;

    // Pick left + right primers (standard PCR pair)
    pa->pick_left_primer  = 1;
    pa->pick_right_primer = 1;
    pa->pick_internal_oligo = 0;

    // --- 2. Sequence arguments ---
    seq_args *sa = create_seq_arg();

    const char *tmpl = ctx.tmpl.c_str();

    sa->sequence = strdup(tmpl);
    sa->sequence_name = strdup("my_template");
    
    sa->ok_regions.left_pairs[0][1]  = ctx.args.len_PDR;
    sa->ok_regions.right_pairs[0][1] = ctx.args.len_PDR;
    sa->ok_regions.any_left  = 0;
    sa->ok_regions.any_right = 0;
    sa->ok_regions.any_pair  = 0;
    sa->ok_regions.count     = 1;


    // sa->incl_s = ctx.pdr_regions[0] + ctx.args.len_PDR;
    // sa->incl_l = ctx.pdr_regions[1] - sa->incl_s;

    // Number of primer pairs to return
    pa->num_return = 20;

    // --- 3. Run primer design ---
    for (std::size_t i = 0; i < ctx.pdr_regions.size(); i += 2) {

        sa->ok_regions.left_pairs[0][0]  = ctx.pdr_regions[i];
        sa->ok_regions.right_pairs[0][0] = ctx.pdr_regions[i + 1];

        p3retval *retval = choose_primers(pa, sa);

        if (retval == nullptr) {
            std::cerr << "choose_primers() returned NULL\n";
            return ;
        }

        // Check for errors
        if (retval->glob_err.data) {
            std::cerr << "Global error: " << retval->glob_err.data << "\n";
        }
        if (retval->per_sequence_err.data) {
            std::cerr << "Sequence error: " << retval->per_sequence_err.data << "\n";
        }

        // --- 4. Extract results ---
        int n_pairs = retval->best_pairs.num_pairs;
        std::cout << "PDR pair " << i / 2 << ": "
                  << "Found " << n_pairs << " primer pair(s)\n";
        ctx.candidate_primers.push_back(extract_all(retval, sa));

        destroy_p3retval(retval);
    }

    destroy_seq_args(sa);
    p3_destroy_global_settings(pa);
}

void DimerStage::run(PipelineContext& ctx) {
    srand(ctx.args.seed);
    /*
    KPartiteGraph g(ctx.args.input_file);
    auto solution = g.solve_fast(ctx.args.iter);
    std::cout << "loss: " << g.cost(solution) << std::endl;

    std::ofstream fout(ctx.args.output_file + "." + short_name());
    for (index_t i : solution)
        fout << i << " ";
    fout << "\n";
    fout.close();

    ctx.dimer_solution = solution;
    */
}

void OffTargetStage::run(PipelineContext& ctx) {
    /*
    std::cout << "Running off-target search\n";
    Thal::init(std::string(PRIMER3_PATH) + "/src/primer3_config", 
            args.mv, 
            args.dv, 
            args.dntp, 
            args.dna_conc, 
            args.temp);

    std::vector<std::string> labels, data;
    read_fasta(args.input_file, labels, data);

    Automaton *ac = new Automaton(labels, data, args.kmer_len);
    auto results = ac->search(args.ref_file, 
            args.kmer_len - 1, 
            args.threshold,
            args.dg_thres,
            args.chunk_size, 
            args.block_size, 
            args.nthreads);
    delete ac;
    
    write_results(args.output_file, results, labels);
    */
}

void Pipeline::run(PipelineContext& ctx) {
    std::cout << "\n" 
              << std::string(50, '=') << "\n"
              << std::string(11, ' ') << "PIPELINE EXECUTION STARTED\n"
              << std::string(50, '=') << "\n";
    
    for (size_t i = 0; i < stages_.size(); ++i) {
        std::cout << "\n[" << (i+1) << "/" << stages_.size() << "] "
                  << "Running: " << stages_[i]->name() << "\n"
                  << std::string(50, '-') << "\n";
        
        stages_[i]->run(ctx);
        
        std::cout << "✓ Completed: " << stages_[i]->name() << "\n"
                  << std::string(50, '-') << "\n";
    }
}
