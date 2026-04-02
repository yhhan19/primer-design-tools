#include "utility.hpp"
#include "automaton.hpp"
#include "risk_optimizer.hpp"
#include "graph.hpp"
#include "pipeline.hpp"

static void print_usage(const char* prog) {
    std::cerr
        << "Usage:\n"
        << "  " << prog << " <mode> -i <input> -o <output> [options]\n"
        << "\n"
        << "Modes (required):\n"
        << "  pdr               Run PDR optimization\n"
        << "  dimer             Run dimer optimization\n"
        << "  off               Run off-target search\n"
        << "\n"
        << "I/O (required):\n"
        << "  -i <file>         Input file\n"
        << "  -o <file>         Output file\n"
        << "  -r <file>         Reference file (required for off mode)\n"
        << "\n"
        << "PDR / risk optimization parameters (mode = pdr):\n"
        << "  -S <int>          Random seed (default 42)\n"
        << "  -Ln <int>         Amplicon max length (default 420)\n"
        << "  -Lx <int>         Amplicon min length (default 252)\n"
        << "  -Lp <int>         PDR length (default 40)\n"
        << "  -Ux <double>      Max risk (default 10000)\n"
        << "  -Un <double>      Min risk (default 0.1)\n"
        << "\n"
        << "Dimer parameters (mode = dimer):\n"
        << "  -I <int>          Iterations for solver (default 1000)\n"
        << "\n"
        << "Off-target search parameters (mode = off):\n"
        << "  -k <int>          K-mer length (default 8)\n"
        << "  -H <int>          Hamming threshold (default 2)\n"
        << "  -G <double>       Delta G threshold (default 1e8)\n"
        << "  -t <int>          Number of threads (default 8)\n"
        << "  -C <int>          Chunk size (default 100000000)\n"
        << "  -B <int>          Block size (default 134217728)\n"
        << "\n"
        << "Thermodynamic parameters (mode = off):\n"
        << "  --mv <double>     Monovalent conc (default 50.0)\n"
        << "  --dv <double>     Divalent conc (default 1.5)\n"
        << "  --dntp <double>   dNTP conc (default 0.6)\n"
        << "  --dna <double>    DNA conc (default 200.0)\n"
        << "  --temp <double>   Temperature in C (default 37.0)\n"
        << "\n"
        << "Other:\n"
        << "  -h, --help        Show this help message\n"
        << "\n"
        << "Examples:\n"
        << "  " << prog << " pdr   -i in.txt -o out.txt\n"
        << "  " << prog << " dimer -i in.txt -o out.txt\n"
        << "  " << prog << " off   -i primers.fa -r ref.fa -o hits.txt -t 16\n"
        << std::endl;
}

static bool is_flag(const char* s, 
                    const char* shortf, 
                    const char* longf = nullptr) {
    if (!s) return false;
    if (shortf && std::strcmp(s, shortf) == 0) return true;
    if (longf && std::strcmp(s, longf) == 0) return true;
    return false;
}

static const char* require_value(int& i, 
                                 int argc, 
                                 char** argv, 
                                 const char* opt) {
    if (i + 1 >= argc) {
        throw std::runtime_error(
            std::string("Missing value for option: ") + opt
        );
    }
    return argv[++i];
}

static std::size_t require_integer(int& i, 
                                   int argc, 
                                   char** argv, 
                                   const char* opt) {
    const char* v = require_value(i, argc, argv, opt);
    try {
        std::size_t pos = 0;
        std::size_t tmp = std::stoull(v, &pos);
        if (pos != std::strlen(v)) {
            throw std::runtime_error("Seed contains non-numeric characters");
        }
        return tmp;
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("Invalid number for ") + opt + ": " + 
            v + " (" + e.what() + ")"
        );
    }
}

static double require_numeric(int& i, 
                              int argc, 
                              char** argv, 
                              const char* opt) {
    const char* v = require_value(i, argc, argv, opt);
    try {
        std::size_t pos = 0;
        double tmp = std::stod(v, &pos);
        if (pos != std::strlen(v)) {
            throw std::runtime_error("contains non-numeric characters");
        }
        return tmp;
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("Invalid number for ") + opt + ": " + 
            v + " (" + e.what() + ")"
        );
    }
}

static Args parse_args(int argc, char** argv) {
    std::cout << "\nCommand: \n";
    for (int i = 0; i < argc; ++i) {
        std::cout << argv[i] << " ";
    }
    std::cout << "\n";

    Args a;

    if (argc < 2) {
        std::cerr << "Error: missing mode (pdr / dimer / off)\n";
        print_usage(argv[0]);
        a.mode = "none";
        return a;
    }

    a.mode = argv[1];

    for (int i = 2; i < argc; ++i) {
        const char* cur = argv[i];

        if (is_flag(cur, "-h", "--help")) {
            a.help = true;
            return a;
        } else if (is_flag(cur, "-i")) {
            a.input_file = require_value(i, argc, argv, "-i");
        } else if (is_flag(cur, "-o")) {
            a.output_file = require_value(i, argc, argv, "-o");
        } else if (is_flag(cur, "-r")) {
            a.ref_file = require_value(i, argc, argv, "-r");
        } else if (is_flag(cur, "-k")) {
            a.kmer_len = require_integer(i, argc, argv, "-k");
        } else if (is_flag(cur, "-H")) {
            a.threshold = require_integer(i, argc, argv, "-H") * 2;
        } else if (is_flag(cur, "-H2")) {
            a.threshold = require_integer(i, argc, argv, "-H2") * 2 + 1;
        } else if (is_flag(cur, "-t")) {
            a.nthreads = require_integer(i, argc, argv, "-t");
        } else if (is_flag(cur, "-C")) {
            a.chunk_size = require_integer(i, argc, argv, "-C");
        } else if (is_flag(cur, "-B")) {
            a.block_size = require_integer(i, argc, argv, "-B");
        } else if (is_flag(cur, "-G")) {
            a.dg_thres = require_numeric(i, argc, argv, "-G");
        } else if (is_flag(cur, "--mv")) {
            a.mv = require_numeric(i, argc, argv, "--mv");
        } else if (is_flag(cur, "--dv")) {
            a.dv = require_numeric(i, argc, argv, "--dv");
        } else if (is_flag(cur, "--dntp")) {
            a.dntp = require_numeric(i, argc, argv, "--dntp");
        } else if (is_flag(cur, "--dna")) {
            a.dna_conc = require_numeric(i, argc, argv, "--dna");
        } else if (is_flag(cur, "--temp")) {
            a.temp = require_numeric(i, argc, argv, "--temp");
        } else if (is_flag(cur, "-S")) {
            a.seed = require_integer(i, argc, argv, "-S");
        } else if (is_flag(cur, "-Ln")) {
            a.len_amp = require_integer(i, argc, argv, "-Ln");
        } else if (is_flag(cur, "-Lx")) {
            a.len_amp_min = require_integer(i, argc, argv, "-Lx");
        } else if (is_flag(cur, "-Lp")) {
            a.len_PDR = require_integer(i, argc, argv, "-Lp");
        } else if (is_flag(cur, "-Ux")) {
            a.u_max = require_numeric(i, argc, argv, "-Ux");
        } else if (is_flag(cur, "-Un")) {
            a.u_min = require_numeric(i, argc, argv, "-Un");
        } else if (is_flag(cur, "-I")) {
            a.iter = require_integer(i, argc, argv, "-I");
        } else {
            throw std::runtime_error(std::string("Unknown option: ") + cur);
        }
    }

    if (a.input_file.empty() || a.output_file.empty()) {
        throw std::runtime_error("Missing required -i / -o");
    }

    if (a.mode == "off" && a.ref_file.empty()) {
        throw std::runtime_error("Missing required -r");
    }

    std::cout << "\nProgram Parameters\n";
    std::cout << "  Mode            : " << a.mode        << "\n";
    std::cout << "  Input file      : " << a.input_file  << "\n";
    std::cout << "  Output file     : " << a.output_file << "\n";

    if (!a.ref_file.empty())
        std::cout << "  Reference file  : " << a.ref_file    << "\n";

    if (a.mode == "pdr") {
        std::cout << "\nPDR Parameters\n";
        std::cout << "  Amp length max  : " << a.len_amp     << "\n";
        std::cout << "  Amp length min  : " << a.len_amp_min << "\n";
        std::cout << "  PDR length      : " << a.len_PDR     << "\n";
        std::cout << "  Max risk        : " << a.u_max       << "\n";
        std::cout << "  Min risk        : " << a.u_min       << "\n";
        std::cout << "  Random seed     : " << a.seed        << "\n";
    }

    else if (a.mode == "dimer") {
        std::cout << "\nDimer Parameters\n";
        std::cout << "  Iterations      : " << a.iter        << "\n";
        std::cout << "  Random seed     : " << a.seed        << "\n";
    }

    else if (a.mode == "off") {
        std::cout << "\nOff-target Search Parameters\n";
        std::cout << "  Kmer length     : " << a.kmer_len    << "\n";
        std::cout << "  Hamming thres   : " << a.threshold/2 << "\n";
        std::cout << "  Threads         : " << a.nthreads    << "\n";
        std::cout << "  Chunk size      : " << a.chunk_size  << "\n";
        std::cout << "  Block size      : " << a.block_size  << "\n";

        std::cout << "\nThermo parameters\n";
        std::cout << "  dG thres        : " << a.dg_thres    << "\n";
        std::cout << "  Monovalent (mv) : " << a.mv          << "\n";
        std::cout << "  Divalent (dv)   : " << a.dv          << "\n";
        std::cout << "  dNTP            : " << a.dntp        << "\n";
        std::cout << "  DNA conc        : " << a.dna_conc    << "\n";
        std::cout << "  Temp (C)        : " << a.temp        << "\n";
    }

    std::cout << std::endl;

    return a;
}

void test_primer3() {
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
    pa->pr_min[0] = 100;
    pa->pr_max[0] = 300;
    pa->num_intervals = 1;

    // Pick left + right primers (standard PCR pair)
    pa->pick_left_primer  = 1;
    pa->pick_right_primer = 1;
    pa->pick_internal_oligo = 0;

    // --- 2. Sequence arguments ---
    seq_args *sa = create_seq_arg();

    const char *tmpl = "GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCC"
                       "TACATTTTAGCATCAGTGAGTACAGCATGCTTACTGG"
                       "AAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTG"
                       "CAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCAT"
                       "GGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTG";

    sa->sequence = strdup(tmpl);
    sa->sequence_name = strdup("my_template");

    // Included region: start=0, length=full sequence
    sa->incl_s = 0;
    sa->incl_l = strlen(tmpl);

    // Number of primer pairs to return
    pa->num_return = 5;

    // --- 3. Run primer design ---
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
    std::cout << "Found " << n_pairs << " primer pair(s)\n\n";

    for (int i = 0; i < n_pairs; i++) {
        primer_pair *pp = &retval->best_pairs.pairs[i];
        primer_rec  *left  = pp->left;
        primer_rec  *right = pp->right;

        // Extract left primer sequence from template
        char left_seq[64] = {0};
        strncpy(left_seq, sa->sequence + left->start, left->length);

        // Right primer is on reverse complement — extract raw coords
        char right_seq[64] = {0};
        // right->start is the 3' end position; length goes leftward
        int right_start = right->start - right->length + 1;
        strncpy(right_seq, sa->sequence + right_start, right->length);

        std::cout << "Pair " << i + 1 << ":\n"
                  << "  Left:  " << left_seq
                  << "  (pos=" << left->start
                  << " len=" << left->length
                  << " Tm=" << left->temp << ")\n"
                  << "  Right: " << right_seq
                  << "  (pos=" << right->start
                  << " len=" << right->length
                  << " Tm=" << right->temp << ")\n"
                  << "  Product size: " << pp->product_size << "\n\n";
    }

    // --- 5. Cleanup ---
    destroy_p3retval(retval);
    destroy_seq_args(sa);
    p3_destroy_global_settings(pa);
}

int main(int argc, char** argv) {
    try {
        PipelineContext ctx;

        ctx.args = parse_args(argc, argv);
        
        if (ctx.args.help) {
            print_usage(argv[0]);
            return 0;
        }

        read_fasta(ctx.args.input_file, ctx.labels, ctx.sequences);
        
        Pipeline p;
        p.add(std::make_unique<PDRStage>());
        p.add(std::make_unique<PrimerSelStage>());
        p.add(std::make_unique<DimerStage>());
        p.add(std::make_unique<OffTargetStage>());
        p.run(ctx);
        
        /*
        if (args.mode == "pdr") {
            std::cout << "Running PDR optimizer\n";
            srand(args.seed);
            std::string ref;
            std::vector<risk_t> rate;
            if (read_rate_ref(args.input_file, rate, ref) != 0) return 1;

            RiskOptimizer ro(rate, 
                             args.len_PDR, 
                             args.len_amp, 
                             args.len_amp);
            auto PDR = ro.search(0, args.u_max, args.u_min);
            ro.validate_PDR(PDR);
            risk_t score = ro.score(PDR, ALPHA);
            std::cout << "loss: " << score << std::endl;

            std::ofstream fout(args.output_file);
            fout << display_PDR(PDR) << std::endl;
        }
        else if (args.mode == "dimer") {
            std::cout << "Running dimer optimization\n";
            srand(args.seed);
            KPartiteGraph g(args.input_file);
            auto solution = g.solve_fast(args.iter);
            std::cout << "loss: " << g.cost(solution) << std::endl;

            std::ofstream fout(args.output_file);
            for (index_t i : solution)
                fout << i << " ";
            fout << "\n";
        }
        else if (args.mode == "off") {
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
        }
        */
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n\n";
        print_usage(argv[0]);
        return 2;
    }
}
