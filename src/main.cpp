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

int main(int argc, char** argv) {
    try {
        PipelineContext ctx;

        ctx.args = parse_args(argc, argv);
        
        if (ctx.args.help) {
            print_usage(argv[0]);
            return 0;
        }

        read_fasta(ctx.args.input_file, ctx.labels, ctx.sequences);
        ctx.tmpl = msa_consensus(ctx.sequences);
        
        Pipeline p;
        p.add(std::make_unique<PDRStage>());
        p.add(std::make_unique<PrimerSelStage>());
        p.add(std::make_unique<DimerStage>());
        p.add(std::make_unique<OffTargetStage>());
        p.run(ctx);
        
        for (auto out : ctx.candidate_primers) 
            display_primer_output(out);
        
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n\n";
        print_usage(argv[0]);
        return 2;
    }
}
