#include "utility.hpp"
#include "automaton.hpp"
#include "risk_optimizer.hpp"
#include "graph.hpp"
#include "pipeline.hpp"
#include "run_log.hpp"

static void print_usage(const char* prog, bool all = false) {
    const Args defaults;
    std::ostringstream usage;
    usage
        << "Usage:\n"
        << "  " << prog << " -i <input> -o <output> [options]\n"
        << "\n"
        << "I/O and configuration:\n"
        << "  -i <file>         Input FASTA file or directory of .fa/.fasta files (required)\n"
        << "  -o <file>         Output file prefix (required)\n"
        << "  -r <file>         Optional background reference; omit to skip off-target search\n"
        << "  -a <file>         Optional per-position risk file; run PDR stage only\n"
        << "  --args <file>     Load saved key = value configuration\n"
        << "  --log-file <file>  Full plain-text run log (default <output>.log)\n"
        << "  --console-mode <auto|live|plain>  Collapse completed blocks in live terminals; auto uses plain text when redirected\n"
        << "                    Options/configurations apply in command-line order; later values override earlier ones.\n"
        << "                    Effective settings are saved to <output>.args.\n"
        << "\n"
        << "PDR / risk optimization parameters:\n"
        << "  --amplicon-max, -Ln <int>         Amplicon max length (default " << defaults.len_amp << ")\n"
        << "  --amplicon-min, -Lx <int>         Amplicon min length (default " << defaults.len_amp_min << ")\n"
        << "  --pdr-length, -Lp <int>         PDR length (default " << defaults.len_PDR << ")\n"
        << "  --pdr-search <fixed|midpoint>  PDR layout algorithm (default " << defaults.pdr_search << ")\n"
        << "  --pdr-min-candidates <int>  Exact Primer3 candidates required per PDR side (default " << defaults.pdr_min_candidates << ": adaptive)\n"
        << "  -Ux <double>      Max PDR score (default " << defaults.u_max << ")\n"
        << "  -Un <double>      Min PDR score (default " << defaults.u_min << ")\n"
        << "\n"
        << "Primer3 parameters:\n"
        << "  --p3-opt-size <int>     Primer optimal size (default " << defaults.p3_opt_size << ")\n"
        << "  --p3-min-size <int>     Primer min size (default " << defaults.p3_min_size << ")\n"
        << "  --p3-max-size <int>     Primer max size (default " << defaults.p3_max_size << ")\n"
        << "  --p3-opt-tm <double>    Primer optimal Tm (default " << defaults.p3_opt_tm << ")\n"
        << "  --p3-min-tm <double>    Primer min Tm (default " << defaults.p3_min_tm << ")\n"
        << "  --p3-max-tm <double>    Primer max Tm (default " << defaults.p3_max_tm << ")\n"
        << "  --p3-min-gc <double>    Primer min GC% (default " << defaults.p3_min_gc << ")\n"
        << "  --p3-max-gc <double>    Primer max GC% (default " << defaults.p3_max_gc << ")\n"
        << "  --candidate-min-return <int>   Shortage warning / single-mode retention floor (default " << defaults.num_return << ")\n"
        << "\nMSA pair filtering (after pair geometry checks):\n"
        << "  --pair-keep-percent <0-100>  Ranked pairs to retain, including ties (default " << defaults.pair_keep_percent << ")\n"
        << "  --pair-min-return <int>     Minimum pairs retained when available (default " << defaults.pair_min_return << ")\n"
        << "  --pair-max-mismatches <int>  Maximum mismatches per primer for joint coverage (default " << defaults.pair_max_mismatches << ")\n"
        << "  --pair-max-3prime-mismatches <int>  Maximum 3-prime mismatches per primer (default " << defaults.pair_max_3prime_mismatches << ")\n"
        << "  --pair-3prime-window <int>  Terminal window length in bases (default " << defaults.pair_3prime_window << ")\n"
        << "  Rank: joint coverage, joint exact 3-prime match, then weighted mismatch.\n"
        << "  Fractions use all MSA sequences; gaps/ambiguous sites are unknown, not matches.\n"
        << "\n"
        << "Off-target search parameters:\n"
        << "  --off-target-mode <single|paired>  Filtering unit (default " << defaults.off_target_mode << ")\n"
        << "  --off-target-amp-min <int>  Minimum predicted off-target product, bp (default " << defaults.off_target_amp_min << ")\n"
        << "  --off-target-amp-max <int>  Maximum predicted off-target product, bp (default " << defaults.off_target_amp_max << ")\n"
        << "  Paired mode checks retained candidate pairs: same contig, inward-facing sites, both dG < -G.\n"
        << "  Chunk overlap is max product length minus one; the left site owns each product.\n"
        << "  Reference products are not automatically excluded as intended targets.\n"
        << "  -k <int>          Maximum exact seed length (default " << defaults.kmer_len << ")\n"
        << "  -H <int>          Hamming threshold (default " << defaults.threshold/2 << ")\n"
        << "  -H2 <int>         Alternative adjacent-mismatch count (opt-in; same limit units as -H)\n"
        << "  -G <double>       Off-target dG cutoff, cal/mol (default " << defaults.dg_thres << ")\n"
        << "  -t <int>          Number of threads (default " << defaults.nthreads << ")\n"
        << "  -C <int>          Chunk size (default " << defaults.chunk_size << ")\n"
        << "  -B <int>          Block size (default " << defaults.block_size << ")\n"
        << "\n"
        << "Dimer parameters:\n"
        << "  --dimer-methods <list>  Comma-separated sum,min-max,top-k-bound, or one cutoff method alone (default " << defaults.dimer_methods << ")\n"
        << "                           Cutoff methods: min-below-cutoff, min-affected-amplicons, min-affected-primers.\n"
        << "                           min-affected-amplicons minimizes affected targets, then total cutoff deficit.\n"
        << "                           min-affected-primers minimizes affected primer count only; no tie-breakers.\n"
        << "  --dimer-select <method>  Method supplying final panel (default first listed method)\n"
        << "  --dimer-rounds <int>  Reserved; refinement is disabled, only 1 is accepted\n"
        << "  --dimer-restarts <int>  Main selection restarts per method/pool (default " << defaults.dimer_restarts << ")\n"
        << "  --dimer-landscape-points <int>  Independent cutoff optimizations, 0 disables; otherwise 2..100 (default 0)\n"
        << "  --dimer-landscape-adaptive  Allocate cutoff searches adaptively around count changes\n"
        << "  --landscape-only  Run cutoff searches without main selection or selected-panel outputs\n"
        << "  --dimer-landscape-restarts <int>  Restarts per landscape cutoff (default 200)\n"
        << "  --dimer-landscape-iterations <int>  Proposals per restart per landscape cutoff (default 1000)\n"
        << "  --dimer-landscape-min <double>  Lowest sampled cutoff in cal/mol (default -5000)\n"
        << "  --dimer-landscape-max <double>  Highest sampled cutoff in cal/mol (default -2000)\n"
        << "  --dimer-start-prefix <path>  Seed from <path>.poolN.<method>.primers.tsv; candidates must match\n"
        << "  --dimer-cache <directory|auto|off>  Reuse interaction energies; auto uses <output-directory>/.dpro-dimer-cache (default auto)\n"
        << "  --dimer-bad-percent <double>  Allowed worst interactions in percent [0,100), default " << defaults.dimer_bad_percent << "\n"
        << "  --dimer-dg-cutoff <double>  Minimize interactions strictly below this dG (cal/mol; default "
        << Args::default_dimer_dg_cutoff << "); errors with incompatible objectives\n"
        << "  --dimer-temp-min <double>  Annealing window minimum in C (default " << defaults.dimer_temp_min.value() << ")\n"
        << "  --dimer-temp-max <double>  Annealing window maximum in C (default " << defaults.dimer_temp_max.value() << ")\n"
        << "                           Use equal min/max values for a single dimer temperature.\n"
        << "  -I <int>          Solver budget; main selection uses up to 5 * -I proposals/restart (default " << defaults.iter << ")\n"
        << "  -S <int>          Random seed (default " << defaults.seed << ")\n"
        << "\n"
        << "Thermodynamic parameters:\n"
        << "  --mv <double>     Monovalent conc, mM (default " << defaults.mv << ")\n"
        << "  --dv <double>     Divalent conc, mM (default " << defaults.dv << ")\n"
        << "  --dntp <double>   dNTP conc, mM (default " << defaults.dntp << ")\n"
        << "  --dna <double>    DNA conc, nM (default " << defaults.dna_conc << ")\n"
        << "  --temp <double>   Non-dimer scoring temperature, C (default " << defaults.temp << ")\n"
        << "                    Dimer min/max override --temp only for dimer scoring.\n"
        << "                    Retained oligos are paired before off-target screening.\n"
        << "                    Midpoint mode filters actual product length; repeated Primer3 pair checks are disabled.\n"
        << "                    All dimer methods select only validated pair combinations.\n"
        << "\nFixed internal settings (reported in the run log; no CLI switches):\n"
        << "  PDR forward/reverse limits: " << defaults.len_fu << "/" << defaults.len_ru << " bp\n"
        << "  Dimer SA: configurable restarts, 5 * -I proposals/restart, initial temperature 1000 (1 for cutoff counts),\n"
        << "            exponential cooling toward 1e-6; each method resets to -S.\n"
        << "  All dimer objectives exclude self-dimers. Higher dG is better.\n"
        << "  top-k-bound uses k=floor(edges * percent/100), then maximizes rank k+1.\n"
        << "  Final dimer labels use --dimer-dg-cutoff, or -3000 cal/mol for other objectives.\n"
        << "  BED columns 8-12 contain dimer status, partner count, worst bad dG, cutoff, partner names.\n"
        << "  .bad-dimers.tsv lists selected pairs below that cutoff; orange primer outlines mark them in the SVG.\n"
        << "  Saved threshold_mode is standard (-H) or adjacent (-H2); threshold stores the user limit.\n"
        << "\nOptional region expansion:\n"
        << "  --expand-affected-pdrs  Expand affected regions, filter new pairs and rerun selection (default on in full midpoint mode)\n"
        << "  --no-expand-affected-pdrs  Disable expansion explicitly\n"
        << "  --pdr-expansion-bases <int>  Extra bases per edge of affected PDRs, bounded by alternating-pool layout (default 20)\n"
        << "  Original candidates are retained; final layout and configured objective are checked.\n"
        << "\nOther:\n"
        << "  -h, --help        Show common options; --help-all shows all options\n"
        << "\n";
    if (all) {
        std::cerr << usage.str() << "  --p3-num-return <int>  Deprecated alias for --candidate-min-return\n";
        return;
    }
    // Keep common controls concise; the full inventory remains available above.
    const std::vector<std::string> advanced = {
        "  -a ", "--pdr-search", "--pdr-min-candidates", "  -Ux ", "  -Un ",
        "--candidate-min-return", "--dimer-select", "--dimer-rounds", "--dimer-start-prefix",
        "--dimer-bad-percent", "  -H2 ", "  -k ", "  -C ", "  -B "
    };
    std::istringstream lines(usage.str());
    std::string line;
    while (std::getline(lines, line)) {
        if (line == "Fixed internal settings (reported in the run log; no CLI switches):") {
            while (std::getline(lines, line) && !line.empty()) {}
            continue;
        }
        bool hidden = false;
        for (const auto& flag : advanced)
            if (line.find(flag) != std::string::npos) hidden = true;
        if (!hidden) std::cerr << line << '\n';
    }
    std::cerr << "\nAdvanced algorithms (fixed layout, single off-target mode and comparison objectives)\n"
              << "and compatibility aliases are documented in --help-all.\n";

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
        if (v[0] == '-') throw std::runtime_error("must be nonnegative");
        std::size_t pos = 0;
        std::size_t tmp = std::stoull(v, &pos);
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

struct ParseState {
    bool expansion_explicit = false;
};

static Args parse_args_file(const std::string& path, Args a, ParseState& state) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open args file: " + path);

    std::string line;
    while (std::getline(f, line)) {
        // strip comments and blank lines
        auto pos = line.find('#');
        if (pos != std::string::npos) line = line.substr(0, pos);
        if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

        std::istringstream ss(line);
        std::string key, eq, val;
        if (!(ss >> key >> eq) || eq != "=")
            throw std::runtime_error("Malformed args file line: " + line);
        std::getline(ss, val);
        const auto first_value = val.find_first_not_of(" \t\r\n");
        if (first_value == std::string::npos)
            throw std::runtime_error("Missing args file value: " + key);
        val = val.substr(first_value, val.find_last_not_of(" \t\r\n") - first_value + 1);

        if      (key == "input_file")   a.input_file   = val;
        else if (key == "output_file")  a.output_file  = val;
        else if (key == "log_file") a.log_file = val;
        else if (key == "console_mode") a.console_mode = val;
        else if (key == "rate_file") a.rate_file = val;
        else if (key == "ref_file")     a.ref_file     = val;
        else if (key == "seed")         a.seed         = std::stoull(val);
        else if (key == "len_amp")      a.len_amp      = std::stoull(val);
        else if (key == "len_amp_min")  a.len_amp_min  = std::stoull(val);
        else if (key == "len_PDR")      a.len_PDR      = std::stoull(val);
        else if (key == "pdr_search")   a.pdr_search   = val;
        else if (key == "pdr_min_candidates") a.pdr_min_candidates = std::stoull(val);
        else if (key == "u_max")        a.u_max        = std::stod(val);
        else if (key == "u_min")        a.u_min        = std::stod(val);
        else if (key == "iter")         a.iter         = std::stoull(val);
        else if (key == "off_target_mode") a.off_target_mode = val;
        else if (key == "off_target_amp_min") {
            std::size_t consumed=0;
            if(val.empty() || val[0]=='-') throw std::runtime_error("Invalid --off-target-amp-min in args file");
            a.off_target_amp_min=std::stoull(val,&consumed);
            if(consumed!=val.size()) throw std::runtime_error("Invalid --off-target-amp-min in args file");
        }
        else if (key == "off_target_amp_max") {
            std::size_t consumed=0;
            if(val.empty() || val[0]=='-') throw std::runtime_error("Invalid --off-target-amp-max in args file");
            a.off_target_amp_max=std::stoull(val,&consumed);
            if(consumed!=val.size()) throw std::runtime_error("Invalid --off-target-amp-max in args file");
        }
        else if (key == "kmer_len")     a.kmer_len     = std::stoull(val);
        else if (key == "threshold")    a.threshold    = std::stoull(val) * 2 + a.threshold % 2;
        else if (key == "threshold_mode") {
            if (val != "standard" && val != "adjacent")
                throw std::runtime_error("threshold_mode must be standard or adjacent");
            a.threshold = (a.threshold / 2) * 2 + (val == "adjacent" ? 1 : 0);
        }
        else if (key == "nthreads")     a.nthreads     = std::stoull(val);
        else if (key == "block_size")   a.block_size   = std::stoull(val);
        else if (key == "chunk_size")   a.chunk_size   = std::stoull(val);
        else if (key == "dg_thres")     a.dg_thres     = std::stod(val);
        else if (key == "mv")           a.mv           = std::stod(val);
        else if (key == "dv")           a.dv           = std::stod(val);
        else if (key == "dntp")         a.dntp         = std::stod(val);
        else if (key == "dna_conc")     a.dna_conc     = std::stod(val);
        else if (key == "temp")         a.temp         = std::stod(val);
        else if (key == "pair_keep_percent") {
            std::size_t consumed = 0;
            if (val.empty() || val[0] == '-')
                throw std::runtime_error("Invalid --pair-keep-percent in args file");
            a.pair_keep_percent = std::stod(val, &consumed);
            if (consumed != val.size())
                throw std::runtime_error("Invalid --pair-keep-percent in args file");
        }
        else if (key == "pair_min_return") {
            std::size_t consumed = 0;
            if (val.empty() || val[0] == '-')
                throw std::runtime_error("Invalid --pair-min-return in args file");
            a.pair_min_return = std::stoull(val, &consumed);
            if (consumed != val.size())
                throw std::runtime_error("Invalid --pair-min-return in args file");
        }
        else if (key == "pair_max_mismatches") {
            std::size_t consumed = 0;
            if (val.empty() || val[0] == '-')
                throw std::runtime_error("Invalid --pair-max-mismatches in args file");
            a.pair_max_mismatches = std::stoull(val, &consumed);
            if (consumed != val.size())
                throw std::runtime_error("Invalid --pair-max-mismatches in args file");
        }
        else if (key == "pair_max_3prime_mismatches") {
            std::size_t consumed = 0;
            if (val.empty() || val[0] == '-')
                throw std::runtime_error("Invalid --pair-max-3prime-mismatches in args file");
            a.pair_max_3prime_mismatches = std::stoull(val, &consumed);
            if (consumed != val.size())
                throw std::runtime_error("Invalid --pair-max-3prime-mismatches in args file");
        }
        else if (key == "pair_3prime_window") {
            std::size_t consumed = 0;
            if (val.empty() || val[0] == '-')
                throw std::runtime_error("Invalid --pair-3prime-window in args file");
            a.pair_3prime_window = std::stoull(val, &consumed);
            if (consumed != val.size())
                throw std::runtime_error("Invalid --pair-3prime-window in args file");
        }
        else if (key == "dimer_methods") a.dimer_methods = val;
        else if (key == "expand_affected_pdrs") {
            if (val != "0" && val != "1") throw std::runtime_error("Invalid expand_affected_pdrs");
            a.expand_affected_pdrs = val == "1";
            state.expansion_explicit = true;
        }
        else if (key == "pdr_expansion_bases") {
            std::size_t parsed = 0;
            if (val.empty() || val[0] == '-') throw std::runtime_error("Invalid pdr_expansion_bases");
            a.pdr_expansion_bases = std::stoull(val, &parsed);
            if (parsed != val.size()) throw std::runtime_error("Invalid pdr_expansion_bases");
        }
        else if (key == "landscape_only") {
            if (val != "0" && val != "1") throw std::runtime_error("Invalid landscape_only");
            a.landscape_only = val == "1";
        }
        else if (key == "dimer_landscape_adaptive") {
            if (val != "0" && val != "1") throw std::runtime_error("Invalid dimer_landscape_adaptive");
            a.dimer_landscape_adaptive = val == "1";
        }
        else if (key == "dimer_start_prefix") a.dimer_start_prefix = val;
        else if (key == "dimer_cache") a.dimer_cache = val;
        else if (key == "dimer_restarts" || key == "dimer_landscape_restarts" || key == "dimer_landscape_iterations") {
            std::size_t parsed = 0;
            if (val.empty() || val[0] == '-') throw std::runtime_error("Invalid " + key);
            const auto value = std::stoull(val, &parsed);
            if (parsed != val.size()) throw std::runtime_error("Invalid " + key);
            if (key == "dimer_restarts") a.dimer_restarts = value;
            else if (key == "dimer_landscape_restarts") a.dimer_landscape_restarts = value;
            else a.dimer_landscape_iterations = value;
        }
        else if (key == "dimer_rounds" || key == "dimer_landscape_points") {
            std::size_t parsed = 0;
            if (val.empty() || val[0] == '-') throw std::runtime_error("Invalid " + key);
            const auto value = std::stoull(val, &parsed);
            if (parsed != val.size()) throw std::runtime_error("Invalid " + key);
            if (key == "dimer_rounds") a.dimer_rounds = value;
            else a.dimer_landscape_points = value;
        }
        else if (key == "dimer_landscape_min" || key == "dimer_landscape_max") {
            std::size_t parsed = 0;
            const double value = std::stod(val, &parsed);
            if (parsed != val.size()) throw std::runtime_error("Invalid " + key);
            if (key == "dimer_landscape_min") a.dimer_landscape_min = value;
            else a.dimer_landscape_max = value;
        }
        else if (key == "dimer_select") a.dimer_select = val;
        else if (key == "dimer_bad_percent") {
            std::size_t parsed = 0;
            a.dimer_bad_percent = std::stod(val, &parsed);
            if (parsed != val.size())
                throw std::runtime_error("Invalid --dimer-bad-percent in args file");
        }
        else if (key == "dimer_dg_cutoff") {
            std::size_t parsed = 0;
            a.dimer_dg_cutoff = std::stod(val, &parsed);
            if (parsed != val.size())
                throw std::runtime_error("Invalid --dimer-dg-cutoff in args file");
        }
        else if (key == "dimer_temp_min") a.dimer_temp_min = std::stod(val);
        else if (key == "dimer_temp_max") a.dimer_temp_max = std::stod(val);
        else if (key == "p3_opt_size")  a.p3_opt_size  = std::stoull(val);
        else if (key == "p3_min_size")  a.p3_min_size  = std::stoull(val);
        else if (key == "p3_max_size")  a.p3_max_size  = std::stoull(val);
        else if (key == "p3_opt_tm")    a.p3_opt_tm    = std::stod(val);
        else if (key == "p3_min_tm")    a.p3_min_tm    = std::stod(val);
        else if (key == "p3_max_tm")    a.p3_max_tm    = std::stod(val);
        else if (key == "p3_min_gc")    a.p3_min_gc    = std::stod(val);
        else if (key == "p3_max_gc")    a.p3_max_gc    = std::stod(val);
        else if (key == "num_return")   a.num_return   = std::stoull(val);
        else if (key == "p3_weighted_keep_percent") {
            std::cerr << "Warning: deprecated config key p3_weighted_keep_percent is ignored\n";
        }
        else throw std::runtime_error("Unknown key in args file: " + key);
    }
    return a;
}

static void save_args(const Args& a, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open args save file: " + path);

    auto now = std::chrono::system_clock::now();
    auto t   = std::chrono::system_clock::to_time_t(now);
    f << "# Saved args -- " << std::ctime(&t);
    f << std::setprecision(std::numeric_limits<double>::max_digits10);

    f << "\n# I/O\n";
    f << "input_file  = " << a.input_file  << "\n";
    f << "output_file = " << a.output_file << "\n";
    if (!a.log_file.empty()) f << "log_file = " << a.log_file << "\n";
    f << "console_mode = " << a.console_mode << "\n";
    if (!a.ref_file.empty())
    f << "ref_file    = " << a.ref_file    << "\n";

    if (!a.rate_file.empty()) f << "rate_file = " << a.rate_file << "\n";

    f << "\n# PDR\n";
    f << "len_amp     = " << a.len_amp     << "\n";
    f << "len_amp_min = " << a.len_amp_min << "\n";
    f << "len_PDR     = " << a.len_PDR     << "\n";
    f << "pdr_search  = " << a.pdr_search  << "\n";
    f << "pdr_min_candidates = " << a.pdr_min_candidates << "\n";
    f << "u_max       = " << a.u_max       << "\n";
    f << "u_min       = " << a.u_min       << "\n";
    f << "seed        = " << a.seed        << "\n";

    f << "\n# MSA pair filtering\n";
    f << "pair_keep_percent = " << a.pair_keep_percent << "\n";
    f << "pair_min_return = " << a.pair_min_return << "\n";
    f << "pair_max_mismatches = " << a.pair_max_mismatches << "\n";
    f << "pair_max_3prime_mismatches = " << a.pair_max_3prime_mismatches << "\n";
    f << "pair_3prime_window = " << a.pair_3prime_window << "\n";

    f << "\n# Dimer\n";
    f << "iter        = " << a.iter        << "\n";
    f << "dimer_restarts = " << a.dimer_restarts << "\n";
    f << "dimer_methods = " << a.dimer_methods << "\n";
    f << "expand_affected_pdrs = " << a.expand_affected_pdrs << "\n";
    f << "pdr_expansion_bases = " << a.pdr_expansion_bases << "\n";
    f << "dimer_rounds = " << a.dimer_rounds << "\n";
    f << "dimer_landscape_points = " << a.dimer_landscape_points << "\n";
    f << "dimer_landscape_adaptive = " << a.dimer_landscape_adaptive << "\n";
    f << "landscape_only = " << a.landscape_only << "\n";
    f << "dimer_landscape_restarts = " << a.dimer_landscape_restarts << "\n";
    f << "dimer_landscape_iterations = " << a.dimer_landscape_iterations << "\n";
    f << std::setprecision(std::numeric_limits<double>::max_digits10);
    f << "dimer_landscape_min = " << a.dimer_landscape_min << "\n";
    f << "dimer_landscape_max = " << a.dimer_landscape_max << "\n";
    if (!a.dimer_start_prefix.empty()) f << "dimer_start_prefix = " << a.dimer_start_prefix << "\n";
    f << "dimer_cache = " << a.dimer_cache << "\n";
    if (!a.dimer_select.empty()) f << "dimer_select = " << a.dimer_select << "\n";
    f << "dimer_bad_percent = " << std::setprecision(std::numeric_limits<double>::max_digits10)
      << a.dimer_bad_percent << "\n";
    if (a.dimer_dg_cutoff) f << "dimer_dg_cutoff = " << *a.dimer_dg_cutoff << "\n";
    if (a.dimer_temp_min) f << "dimer_temp_min = " << *a.dimer_temp_min << "\n";
    if (a.dimer_temp_max) f << "dimer_temp_max = " << *a.dimer_temp_max << "\n";

    f << "\n# Off-target search\n";
    f << "off_target_mode = " << a.off_target_mode << "\n";
    f << "off_target_amp_min = " << a.off_target_amp_min << "\n";
    f << "off_target_amp_max = " << a.off_target_amp_max << "\n";
    f << "kmer_len    = " << a.kmer_len    << "\n";
    f << "threshold   = " << a.threshold/2 << "\n";
    f << "threshold_mode = " << (a.threshold % 2 ? "adjacent" : "standard") << "\n";
    f << "nthreads    = " << a.nthreads    << "\n";
    f << "chunk_size  = " << a.chunk_size  << "\n";
    f << "block_size  = " << a.block_size  << "\n";

    f << "\n# Thermodynamics\n";
    f << "dg_thres    = " << a.dg_thres    << "\n";
    f << "mv          = " << a.mv          << "\n";
    f << "dv          = " << a.dv          << "\n";
    f << "dntp        = " << a.dntp        << "\n";
    f << "dna_conc    = " << a.dna_conc    << "\n";
    f << "temp        = " << a.temp        << "\n";

    f << "\n# Primer3\n";
    f << "p3_opt_size = " << a.p3_opt_size << "\n";
    f << "p3_min_size = " << a.p3_min_size << "\n";
    f << "p3_max_size = " << a.p3_max_size << "\n";
    f << "p3_opt_tm   = " << a.p3_opt_tm   << "\n";
    f << "p3_min_tm   = " << a.p3_min_tm   << "\n";
    f << "p3_max_tm   = " << a.p3_max_tm   << "\n";
    f << "p3_min_gc   = " << a.p3_min_gc   << "\n";
    f << "p3_max_gc   = " << a.p3_max_gc   << "\n";
    f << "num_return  = " << a.num_return  << "\n";
}

static Args parse_args(int argc, char** argv) {
    Args a;
    std::vector<std::string> loaded_args_files;
    ParseState state;
 
    // Print command line start
    std::cout << "\nCommand executed:" << std::endl;
    std::cout << argv[0];
 
    for (int i = 1; i < argc; ++i) {
        const char* cur = argv[i];
 
        // Print current argument with formatting
        std::cout << " \\" << std::endl << "  " << cur;
 
        if (is_flag(cur, "-h", "--help") || is_flag(cur, "--help-all")) {
            a.help_all = is_flag(cur, "--help-all");
            a.help = true;
            std::cout << std::endl << std::endl;
            return a;
        } else if (is_flag(cur, "--args")) {
            const char* path = require_value(i, argc, argv, "--args");
            a = parse_args_file(path, a, state);
            loaded_args_files.push_back(path);
            std::cout << " " << path << "  # Args file";
        } else if (is_flag(cur, "--log-file")) {
            a.log_file = require_value(i, argc, argv, "--log-file");
        } else if (is_flag(cur, "--console-mode")) {
            a.console_mode = require_value(i, argc, argv, "--console-mode");
        } else if (is_flag(cur, "-i")) {
            a.input_file = require_value(i, argc, argv, "-i");
            std::cout << " " << a.input_file << "  # Input file";
        } else if (is_flag(cur, "-o")) {
            a.output_file = require_value(i, argc, argv, "-o");
            std::cout << " " << a.output_file << "  # Output file";
        } else if (is_flag(cur, "-r")) {
            a.ref_file = require_value(i, argc, argv, "-r");
            std::cout << " " << a.ref_file << "  # Reference file";
        } else if (is_flag(cur, "-a")) {
            a.rate_file = require_value(i, argc, argv, "-a");
            std::cout << " " << a.rate_file << "  # Rate file";
        } else if (is_flag(cur, "--off-target-mode")) {
            a.off_target_mode = require_value(i, argc, argv, "--off-target-mode");
            std::cout << " " << a.off_target_mode;
        } else if (is_flag(cur, "--off-target-amp-min")) {
            a.off_target_amp_min = require_integer(i, argc, argv, "--off-target-amp-min");
            std::cout << " " << a.off_target_amp_min;
        } else if (is_flag(cur, "--off-target-amp-max")) {
            a.off_target_amp_max = require_integer(i, argc, argv, "--off-target-amp-max");
            std::cout << " " << a.off_target_amp_max;
        } else if (is_flag(cur, "-k")) {
            a.kmer_len = require_integer(i, argc, argv, "-k");
            std::cout << " " << a.kmer_len << "  # Kmer length";
        } else if (is_flag(cur, "-H")) {
            a.threshold = require_integer(i, argc, argv, "-H") * 2;
            std::cout << " " << (a.threshold/2) << "  # Hamming threshold";
        } else if (is_flag(cur, "-H2")) {
            a.threshold = require_integer(i, argc, argv, "-H2") * 2 + 1;
            std::cout << " " << ((a.threshold-1)/2) << "  # Hamming threshold (odd)";
        } else if (is_flag(cur, "-t")) {
            a.nthreads = require_integer(i, argc, argv, "-t");
            std::cout << " " << a.nthreads << "  # Threads";
        } else if (is_flag(cur, "-C")) {
            a.chunk_size = require_integer(i, argc, argv, "-C");
            std::cout << " " << a.chunk_size << "  # Chunk size";
        } else if (is_flag(cur, "-B")) {
            a.block_size = require_integer(i, argc, argv, "-B");
            std::cout << " " << a.block_size << "  # Block size";
        } else if (is_flag(cur, "-G")) {
            a.dg_thres = require_numeric(i, argc, argv, "-G");
            std::cout << " " << a.dg_thres << "  # dG threshold";
        } else if (is_flag(cur, "--mv")) {
            a.mv = require_numeric(i, argc, argv, "--mv");
            std::cout << " " << a.mv << "  # Monovalent conc, mMentration";
        } else if (is_flag(cur, "--dv")) {
            a.dv = require_numeric(i, argc, argv, "--dv");
            std::cout << " " << a.dv << "  # Divalent conc, mMentration";
        } else if (is_flag(cur, "--dntp")) {
            a.dntp = require_numeric(i, argc, argv, "--dntp");
            std::cout << " " << a.dntp << "  # dNTP conc, mMentration";
        } else if (is_flag(cur, "--dna")) {
            a.dna_conc = require_numeric(i, argc, argv, "--dna");
            std::cout << " " << a.dna_conc << "  # DNA concentration";
        } else if (is_flag(cur, "--pair-keep-percent")) {
            a.pair_keep_percent = require_numeric(i, argc, argv, "--pair-keep-percent");
            std::cout << " " << a.pair_keep_percent;
        } else if (is_flag(cur, "--pair-min-return")) {
            a.pair_min_return = require_integer(i, argc, argv, "--pair-min-return");
            std::cout << " " << a.pair_min_return;
        } else if (is_flag(cur, "--pair-max-mismatches")) {
            a.pair_max_mismatches = require_integer(i, argc, argv, "--pair-max-mismatches");
            std::cout << " " << a.pair_max_mismatches;
        } else if (is_flag(cur, "--pair-max-3prime-mismatches")) {
            a.pair_max_3prime_mismatches = require_integer(i, argc, argv, "--pair-max-3prime-mismatches");
            std::cout << " " << a.pair_max_3prime_mismatches;
        } else if (is_flag(cur, "--pair-3prime-window")) {
            a.pair_3prime_window = require_integer(i, argc, argv, "--pair-3prime-window");
            std::cout << " " << a.pair_3prime_window;
        } else if (is_flag(cur, "--dimer-methods")) {
            a.dimer_methods = require_value(i, argc, argv, "--dimer-methods");
            std::cout << " " << a.dimer_methods;
        } else if (is_flag(cur, "--expand-affected-pdrs")) {
            a.expand_affected_pdrs = true;
            state.expansion_explicit = true;
        } else if (is_flag(cur, "--no-expand-affected-pdrs")) {
            a.expand_affected_pdrs = false;
            state.expansion_explicit = true;
        } else if (is_flag(cur, "--pdr-expansion-bases")) {
            a.pdr_expansion_bases = require_integer(i, argc, argv, "--pdr-expansion-bases");
        } else if (is_flag(cur, "--dimer-rounds")) {
            a.dimer_rounds = require_integer(i, argc, argv, "--dimer-rounds");
        } else if (is_flag(cur, "--dimer-restarts")) {
            a.dimer_restarts = require_integer(i, argc, argv, "--dimer-restarts");
        } else if (is_flag(cur, "--dimer-landscape-points")) {
            a.dimer_landscape_points = require_integer(i, argc, argv, "--dimer-landscape-points");
        } else if (is_flag(cur, "--dimer-landscape-adaptive")) {
            a.dimer_landscape_adaptive = true;
        } else if (is_flag(cur, "--landscape-only")) {
            a.landscape_only = true;
        } else if (is_flag(cur, "--dimer-landscape-restarts")) {
            a.dimer_landscape_restarts = require_integer(i, argc, argv, "--dimer-landscape-restarts");
        } else if (is_flag(cur, "--dimer-landscape-iterations")) {
            a.dimer_landscape_iterations = require_integer(i, argc, argv, "--dimer-landscape-iterations");
        } else if (is_flag(cur, "--dimer-landscape-min")) {
            a.dimer_landscape_min = require_numeric(i, argc, argv, "--dimer-landscape-min");
        } else if (is_flag(cur, "--dimer-landscape-max")) {
            a.dimer_landscape_max = require_numeric(i, argc, argv, "--dimer-landscape-max");
        } else if (is_flag(cur, "--dimer-start-prefix")) {
            a.dimer_start_prefix = require_value(i, argc, argv, "--dimer-start-prefix");
        } else if (is_flag(cur, "--dimer-cache")) {
            a.dimer_cache = require_value(i, argc, argv, "--dimer-cache");
        } else if (is_flag(cur, "--dimer-select")) {
            a.dimer_select = require_value(i, argc, argv, "--dimer-select");
            std::cout << " " << a.dimer_select;
        } else if (is_flag(cur, "--dimer-bad-percent")) {
            a.dimer_bad_percent = require_numeric(i, argc, argv, "--dimer-bad-percent");
            std::cout << " " << a.dimer_bad_percent;
        } else if (is_flag(cur, "--dimer-dg-cutoff")) {
            a.dimer_dg_cutoff = require_numeric(i, argc, argv, "--dimer-dg-cutoff");
            std::cout << " " << *a.dimer_dg_cutoff;
        } else if (is_flag(cur, "--dimer-temp-min")) {
            a.dimer_temp_min = require_numeric(i, argc, argv, "--dimer-temp-min");
            std::cout << " " << *a.dimer_temp_min;
        } else if (is_flag(cur, "--dimer-temp-max")) {
            a.dimer_temp_max = require_numeric(i, argc, argv, "--dimer-temp-max");
            std::cout << " " << *a.dimer_temp_max;
        } else if (is_flag(cur, "--temp")) {
            a.temp = require_numeric(i, argc, argv, "--temp");
            std::cout << " " << a.temp << "  # Temperature";
        } else if (is_flag(cur, "-S")) {
            a.seed = require_integer(i, argc, argv, "-S");
            std::cout << " " << a.seed << "  # Random seed";
        } else if (is_flag(cur, "-Ln", "--amplicon-max")) {
            a.len_amp = require_integer(i, argc, argv, "-Ln");
            std::cout << " " << a.len_amp << "  # Amplicon length max";
        } else if (is_flag(cur, "-Lx", "--amplicon-min")) {
            a.len_amp_min = require_integer(i, argc, argv, "-Lx");
            std::cout << " " << a.len_amp_min << "  # Amplicon length min";
        } else if (is_flag(cur, "-Lp", "--pdr-length")) {
            a.len_PDR = require_integer(i, argc, argv, "-Lp");
            std::cout << " " << a.len_PDR << "  # PDR length";
        } else if (is_flag(cur, "--pdr-search")) {
            a.pdr_search = require_value(i, argc, argv, "--pdr-search");
            std::cout << " " << a.pdr_search << "  # PDR search algorithm";
        } else if (is_flag(cur, "--pdr-min-candidates")) {
            a.pdr_min_candidates = require_integer(
                i, argc, argv, "--pdr-min-candidates");
            std::cout << " " << a.pdr_min_candidates
                      << "  # PDR candidate minimum (0=adaptive)";
        } else if (is_flag(cur, "-Ux")) {
            a.u_max = require_numeric(i, argc, argv, "-Ux");
            std::cout << " " << a.u_max << "  # Max PDR score";
        } else if (is_flag(cur, "-Un")) {
            a.u_min = require_numeric(i, argc, argv, "-Un");
            std::cout << " " << a.u_min << "  # Min PDR score";
        } else if (is_flag(cur, "-I")) {
            a.iter = require_integer(i, argc, argv, "-I");
            std::cout << " " << a.iter << "  # Iterations";
        }  else if (is_flag(cur, "--p3-opt-size")) {
            a.p3_opt_size = require_integer(i, argc, argv, "--p3-opt-size");
            std::cout << " " << a.p3_opt_size << "  # Primer optimal size";
        } else if (is_flag(cur, "--p3-min-size")) {
            a.p3_min_size = require_integer(i, argc, argv, "--p3-min-size");
            std::cout << " " << a.p3_min_size << "  # Primer min size";
        } else if (is_flag(cur, "--p3-max-size")) {
            a.p3_max_size = require_integer(i, argc, argv, "--p3-max-size");
            std::cout << " " << a.p3_max_size << "  # Primer max size";
        } else if (is_flag(cur, "--p3-opt-tm")) {
            a.p3_opt_tm = require_numeric(i, argc, argv, "--p3-opt-tm");
            std::cout << " " << a.p3_opt_tm << "  # Primer optimal Tm";
        } else if (is_flag(cur, "--p3-min-tm")) {
            a.p3_min_tm = require_numeric(i, argc, argv, "--p3-min-tm");
            std::cout << " " << a.p3_min_tm << "  # Primer min Tm";
        } else if (is_flag(cur, "--p3-max-tm")) {
            a.p3_max_tm = require_numeric(i, argc, argv, "--p3-max-tm");
            std::cout << " " << a.p3_max_tm << "  # Primer max Tm";
        } else if (is_flag(cur, "--p3-min-gc")) {
            a.p3_min_gc = require_numeric(i, argc, argv, "--p3-min-gc");
            std::cout << " " << a.p3_min_gc << "  # Primer min GC%";
        } else if (is_flag(cur, "--p3-max-gc")) {
            a.p3_max_gc = require_numeric(i, argc, argv, "--p3-max-gc");
            std::cout << " " << a.p3_max_gc << "  # Primer max GC%";
        } else if (is_flag(cur, "--candidate-min-return", "--p3-num-return")) {
            if (is_flag(cur, "--p3-num-return"))
                std::cerr << "Warning: --p3-num-return is deprecated; use --candidate-min-return\n";
            a.num_return = require_integer(i, argc, argv, "--candidate-min-return");
            std::cout << " " << a.num_return << "  # Minimum after filtering";
        } else {
            throw std::runtime_error(std::string("Unknown option: ") + cur);
        }
    }
 
    // End command line display
    std::cout << std::endl << std::endl;

    if (a.off_target_mode != "single" && a.off_target_mode != "paired")
        throw std::runtime_error("--off-target-mode must be single or paired");
    if (a.off_target_amp_min == 0 || a.off_target_amp_min > a.off_target_amp_max)
        throw std::runtime_error("--off-target-amp-min/max require 1 <= min <= max");
    if (a.chunk_size == 0 || a.block_size == 0 || a.nthreads == 0 || a.kmer_len == 0)
        throw std::runtime_error("-C, -B, -t and -k must be positive");
    if (a.off_target_amp_max - 1 > std::numeric_limits<std::size_t>::max() - a.chunk_size)
        throw std::runtime_error("--off-target-amp-max is too large for chunk overlap");
    if (!std::isfinite(a.dg_thres)) throw std::runtime_error("-G must be finite");
    if (!std::isfinite(a.pair_keep_percent) || a.pair_keep_percent < 0 || a.pair_keep_percent > 100)
        throw std::runtime_error("--pair-keep-percent must be finite and in [0,100]");
    if (a.pair_min_return == 0 || a.pair_3prime_window == 0)
        throw std::runtime_error("--pair-min-return and --pair-3prime-window must be positive");
    if (a.pair_max_3prime_mismatches > a.pair_3prime_window)
        throw std::runtime_error("--pair-max-3prime-mismatches cannot exceed --pair-3prime-window");
    const auto is_cutoff_method = [](const std::string& method) {
        return method == "min-below-cutoff" || method == "min-affected-amplicons" ||
               method == "min-affected-primers";
    };
    if (a.dimer_dg_cutoff && !std::isfinite(*a.dimer_dg_cutoff))
        throw std::runtime_error("--dimer-dg-cutoff must be finite");
    std::set<std::string> seen_methods;
    std::istringstream method_stream(a.dimer_methods);
    std::string method;
    if (a.dimer_methods.empty() || a.dimer_methods.back() == ',')
        throw std::runtime_error("--dimer-methods must be a nonempty comma-separated list");
    while (std::getline(method_stream, method, ',')) {
        if (method != "sum" && method != "min-max" && method != "top-k-bound" &&
            !is_cutoff_method(method))
            throw std::runtime_error("Unknown --dimer-methods entry: " + method);
        if (!seen_methods.insert(method).second)
            throw std::runtime_error("Duplicate --dimer-methods entry: " + method);
    }
    const bool cutoff_objective = is_cutoff_method(a.dimer_methods);
    for (const auto& entry : seen_methods)
        if (is_cutoff_method(entry) && seen_methods.size() != 1)
            throw std::runtime_error("--dimer-methods cutoff objectives must be requested alone");
    if (a.dimer_dg_cutoff && !cutoff_objective)
        throw std::runtime_error("--dimer-dg-cutoff is incompatible with --dimer-methods " + a.dimer_methods);
    if (cutoff_objective && !a.dimer_dg_cutoff)
        a.dimer_dg_cutoff = Args::default_dimer_dg_cutoff;
    if (a.dimer_select.empty())
        a.dimer_select = a.dimer_methods.substr(0, a.dimer_methods.find(','));
    if (!seen_methods.count(a.dimer_select))
        throw std::runtime_error("--dimer-select must be included in --dimer-methods");
    if (a.iter == 0 || a.iter > std::numeric_limits<std::size_t>::max() / 5)
        throw std::runtime_error("-I must be positive and allow multiplication by 5");
    if (!a.dimer_restarts)
        throw std::runtime_error("--dimer-restarts must be positive");
    if (!a.pdr_expansion_bases || a.pdr_expansion_bases > std::size_t(std::numeric_limits<int>::max()))
        throw std::runtime_error("--pdr-expansion-bases must be positive and fit Primer3 coordinates");
    const bool expansion_supported = a.pdr_search == "midpoint" &&
        a.rate_file.empty() && !a.landscape_only;
    if (!state.expansion_explicit)
        a.expand_affected_pdrs = expansion_supported;
    if (a.expand_affected_pdrs && !expansion_supported)
        throw std::runtime_error("--expand-affected-pdrs requires midpoint full-pipeline mode without --landscape-only");
    if (a.dimer_rounds != 1)
        throw std::runtime_error("--dimer-rounds is disabled for now; use 1");
    if (a.dimer_landscape_points == 1 || a.dimer_landscape_points > 100)
        throw std::runtime_error("--dimer-landscape-points must be 0 or 2..100");
    if (a.dimer_landscape_adaptive && !a.dimer_landscape_points)
        throw std::runtime_error("--dimer-landscape-adaptive requires --dimer-landscape-points");
    if (a.landscape_only && (!a.dimer_landscape_points || a.expand_affected_pdrs ||
                            !a.rate_file.empty() || !a.dimer_start_prefix.empty()))
        throw std::runtime_error("--landscape-only requires landscape points, full-pipeline input, expansion disabled and no main-panel seed prefix");
    if (!a.dimer_landscape_restarts || !a.dimer_landscape_iterations)
        throw std::runtime_error("--dimer-landscape-restarts and --dimer-landscape-iterations must be positive");
    if (!std::isfinite(a.dimer_landscape_min) || !std::isfinite(a.dimer_landscape_max) ||
        a.dimer_landscape_min >= a.dimer_landscape_max ||
        !std::isfinite(a.dimer_landscape_max - a.dimer_landscape_min))
        throw std::runtime_error("Landscape requires finite min < max");
    if (!std::isfinite(a.dimer_bad_percent) || a.dimer_bad_percent < 0 || a.dimer_bad_percent >= 100) {
        throw std::runtime_error("--dimer-bad-percent must be finite and in [0,100)");
    }
    if (a.dimer_temp_min.has_value() != a.dimer_temp_max.has_value()) {
        throw std::runtime_error("Specify both --dimer-temp-min and --dimer-temp-max");
    }
    Thal::validate_dimer_range(a.dimer_temp_min.value_or(a.temp),
                              a.dimer_temp_max.value_or(a.temp));
    if (a.input_file.empty() || a.output_file.empty()) {
        throw std::runtime_error("Missing required -i / -o");
    }
    if (a.console_mode != "auto" && a.console_mode != "live" && a.console_mode != "plain")
        throw std::runtime_error("--console-mode must be auto, live or plain");
    if (a.num_return == 0) {
        throw std::runtime_error("--candidate-min-return must be at least 1");
    }
    if (a.pdr_search != "fixed" && a.pdr_search != "midpoint") {
        throw std::runtime_error(
            "--pdr-search must be either fixed or midpoint");
    }
    if (a.len_amp_min > a.len_amp) {
        throw std::runtime_error(
            "amplicon minimum length cannot exceed its maximum length");
    }
    if (a.len_amp_min < a.len_PDR || a.len_amp < a.len_PDR) {
        throw std::runtime_error(
            "amplicon lengths must be at least the PDR length");
    }
    configure_run_log(a.log_file.empty() ? a.output_file + ".log" : a.log_file, a.console_mode);
    std::cout << "\nProgram Parameters\n";
    std::cout << "  Full log        : " << (a.log_file.empty() ? a.output_file + ".log" : a.log_file) << "\n";
    std::cout << "  Console mode    : " << a.console_mode << "\n";
    std::cout << "  Input file      : " << a.input_file  << "\n";
    std::cout << "  Output file     : " << a.output_file << "\n";

    std::cout << "  Reference file  : " << (a.ref_file.empty() ? "(not set; off-target search skipped)" : a.ref_file) << "\n";
    std::cout << "  Rate file (-a)  : " << (a.rate_file.empty() ? "(not set)" : a.rate_file) << "\n";
    std::cout << "  Pipeline mode   : " << (a.landscape_only ? "landscape only" : a.rate_file.empty() ? "full pipeline" : "PDR only") << "\n";
    std::cout << "  Args files      : ";
    if (loaded_args_files.empty()) std::cout << "(none)";
    for (std::size_t i = 0; i < loaded_args_files.size(); ++i)
        std::cout << (i ? ", " : "") << loaded_args_files[i];
    std::cout << "\n";

    std::cout << "\nPDR Parameters\n";
    std::cout << "  Amp length max  : " << a.len_amp     << "\n";
    std::cout << "  Amp length min  : " << a.len_amp_min << "\n";
    std::cout << "  PDR length      : " << a.len_PDR     << "\n";
    std::cout << "  PDR search      : " << a.pdr_search  << "\n";
    std::cout << "  PDR candidate min: ";
    if (a.pdr_min_candidates == 0) std::cout << "adaptive\n";
    else std::cout << a.pdr_min_candidates << "\n";
    std::cout << "  Forward limit   : " << a.len_fu << " bp (fixed internal setting)\n";
    std::cout << "  Reverse limit   : " << a.len_ru << " bp (fixed internal setting)\n";
    std::cout << "  Max PDR score   : " << a.u_max       << "\n";
    std::cout << "  Min PDR score   : " << a.u_min       << "\n";
    std::cout << "  Random seed     : " << a.seed        << "\n";

    std::cout << "\nPrimer3 Parameters\n";
    std::cout << "  Pairing         : retained left x right; geometry + MSA filtering (no repeated Primer3 pair checks)\n";
    std::cout << "  Pair geometry   : "
              << (a.pdr_search == "midpoint" ? "filter actual product length using -Lx/-Ln"
                                             : "prefilter disabled; -Lx/-Ln control PDR layout only")
              << "\n";
    std::cout << "  Pair selection  : Primer3-valid, MSA-retained combinations; whole-pair SA moves\n";
    std::cout << "  Opt size        : " << a.p3_opt_size << "\n";
    std::cout << "  Min size        : " << a.p3_min_size << "\n";
    std::cout << "  Max size        : " << a.p3_max_size << "\n";
    std::cout << "  Opt Tm          : " << a.p3_opt_tm   << "\n";
    std::cout << "  Min Tm          : " << a.p3_min_tm   << "\n";
    std::cout << "  Max Tm          : " << a.p3_max_tm   << "\n";
    std::cout << "  Min GC%         : " << a.p3_min_gc   << "\n";
    std::cout << "  Max GC%         : " << a.p3_max_gc   << "\n";
    std::cout << "  Candidate min   : " << a.num_return  << "\n";
    std::cout << "  Single filtering: disabled; all Primer3 candidates reach pair checks\n";

    std::cout << "\nMSA Pair Filter Parameters\n";
    std::cout << "  Pair keep %     : " << a.pair_keep_percent << " (boundary ties retained)\n";
    std::cout << "  Pair minimum    : " << a.pair_min_return << " (when available)\n";
    std::cout << "  Max mismatches  : " << a.pair_max_mismatches << " per primer\n";
    std::cout << "  Max 3-prime mm  : " << a.pair_max_3prime_mismatches << " per primer\n";
    std::cout << "  3-prime window  : " << a.pair_3prime_window << " bases\n";
    std::cout << "  Pair MSA ranking: joint coverage desc, joint 3-prime match desc, weighted mismatch asc\n";
    std::cout << "  Unknown sites   : gaps/ambiguities; denominator includes all MSA sequences\n";
    std::cout << "\nOff-target Search Parameters\n";
    std::cout << "  Filter mode     : " << a.off_target_mode << "\n";
    std::cout << "  Off-target bp   : " << a.off_target_amp_min << " to " << a.off_target_amp_max << "\n";
    std::cout << "  Paired overlap  : " << a.off_target_amp_max - 1 << " bases\n";
    std::cout << "  Pair criteria   : both dG < cutoff, same contig, inward-facing non-overlapping sites\n";
    std::cout << "  Reference scope : all products; intended loci are not automatically excluded\n";
    std::cout << "  Kmer length     : " << a.kmer_len    << "\n";
    std::cout << "  Hamming thres   : " << a.threshold/2 << "\n";
    std::cout << "  Mismatch mode   : " << (a.threshold % 2 ? "adjacent (-H2)" : "standard (-H)") << "\n";
    std::cout << "  Encoded threshold: " << a.threshold << "\n";
    std::cout << "  Threads         : " << a.nthreads    << "\n";
    std::cout << "  Chunk size      : " << a.chunk_size  << "\n";
    std::cout << "  Block size      : " << a.block_size  << "\n";

    std::cout << "\nDimer Parameters\n";
    std::cout << "  Dimer methods   : " << a.dimer_methods << "\n";
    std::cout << "  Dimer rounds    : " << a.dimer_rounds << "\n";
    std::cout << "  Expand PDRs     : " << a.expand_affected_pdrs << " (+/-"
              << a.pdr_expansion_bases << " bp)\n";
    std::cout << "  Adaptive landscape: " << a.dimer_landscape_adaptive << "\n";
    std::cout << "  Landscape budget: " << a.dimer_landscape_restarts << " restarts x "
              << a.dimer_landscape_iterations << " proposals per cutoff\n";
    std::cout << "  Dimer landscape : " << a.dimer_landscape_points << " points ["
              << a.dimer_landscape_min << ", " << a.dimer_landscape_max << "] cal/mol\n";
    std::cout << "  Dimer start     : " << a.dimer_start_prefix << "\n";
    std::cout << "  Dimer select    : " << (a.dimer_select.empty() ? a.dimer_methods.substr(0, a.dimer_methods.find(',')) : a.dimer_select) << "\n";
    std::cout << "  Dimer bad %     : " << a.dimer_bad_percent << "\n";
    if (a.dimer_dg_cutoff)
        std::cout << "  Dimer cutoff    : " << *a.dimer_dg_cutoff << " cal/mol\n";
    std::cout << "  dG window (C)   : " << a.dimer_temp_min.value_or(a.temp)
              << " to " << a.dimer_temp_max.value_or(a.temp)
              << " (minimum dG, cal/mol)\n";
    std::cout << "  Iterations (-I) : " << a.iter << "\n";
    std::cout << "  SA restarts     : " << a.dimer_restarts << " per method/pool\n";
    std::cout << "  SA proposals    : " << a.iter * 5 << " per restart (5 * -I)\n";
    std::cout << "  SA initial temp : 1000 (score units; fixed)\n";
    std::cout << "  SA final target : 1e-6 (fixed)\n";
    std::cout << "  SA cooling      : " << std::exp(std::log(1e-6 / 1000.0) / double(a.iter * 5)) << " per proposal\n";
    std::cout << "  Dimer rank      : floor(edges * bad_percent / 100) + 1, per pool\n";
    std::cout << "  Self-dimers     : excluded from objectives\n";
    std::cout << "  Random seed     : " << a.seed        << "\n";

    std::cout << "\nThermo parameters\n";
    std::cout << "  dG thres (cal/mol): " << a.dg_thres    << "\n";
    std::cout << "  Monovalent (mM) : " << a.mv          << "\n";
    std::cout << "  Divalent (mM)   : " << a.dv          << "\n";
    std::cout << "  dNTP (mM)       : " << a.dntp        << "\n";
    std::cout << "  DNA conc (nM)   : " << a.dna_conc    << "\n";
    std::cout << "  Temp (C)        : " << a.temp        << "\n";

    std::cout << std::endl;

    std::string save_path = a.output_file + ".args";
    save_args(a, save_path);
    std::cout << "Args saved to: " << save_path << "\n\n";
    return a;
}

int main(int argc, char** argv) {
    bool help_requested = false;
    for (int i=1; i<argc; ++i)
        help_requested |= std::string(argv[i])=="--help-all" || std::string(argv[i])=="--help" || std::string(argv[i])=="-h";
    RunLog run_log(!help_requested);
    try {
        PipelineContext ctx;
        {
            RunLogBlock configuration("Configuration");
            ctx.args = parse_args(argc, argv);
            configuration.complete(ctx.args.help ? "help" : "settings saved");
        }
        
        if (ctx.args.help) {
            print_usage(argv[0], ctx.args.help_all);
            return 0;
        }

        {
        RunLogBlock inputs("Load alignments");
        std::vector<std::string> fasta_files;

        if (fs::is_directory(ctx.args.input_file)) {
            for (auto& entry : fs::directory_iterator(ctx.args.input_file)) {
                auto ext = entry.path().extension().string();
                if (ext == ".fasta" || ext == ".fa") {
                    fasta_files.push_back(entry.path().string());
                }
            }
            std::sort(fasta_files.begin(), fasta_files.end());
        } else {
            // backward-compatible: single FASTA file
            fasta_files.push_back(ctx.args.input_file);
        }

        ctx.segments.resize(fasta_files.size());

        for (std::size_t i = 0; i < fasta_files.size(); ++i) {
            auto& seg = ctx.segments[i];
            seg.name = fs::path(fasta_files[i]).stem().string();

            std::size_t n = read_fasta(fasta_files[i], seg.labels, seg.sequences);
            if (n == 0) {
                throw std::runtime_error("No sequences loaded from: " + fasta_files[i]);
            }

            seg.tmpl = msa_consensus(seg.sequences);
            write_fasta(ctx.args.output_file + "." + seg.name + ".tmpl.fasta",
                        {"template"}, {seg.tmpl}, 80);

            std::cout << "  [" << (i+1) << "/" << fasta_files.size() << "] "
                    << seg.name << " — " << n << (n == 1 ? " sequence, " : " sequences, ")
                    << seg.tmpl.size() << " bp\n";
        }
        inputs.complete(std::to_string(ctx.segments.size()) + " target(s)");
        }

        Pipeline p;

        if (! ctx.args.rate_file.empty()) {
            p.add(std::make_unique<PDRStage>());
            p.run(ctx);
            return 0;
        }

        p.add(std::make_unique<PDRStage>());
        p.add(std::make_unique<PrimerSelStage>());
        p.add(std::make_unique<OffTargetStage>());
        p.add(std::make_unique<DimerStage>());
        if (ctx.args.expand_affected_pdrs) p.add(std::make_unique<PDRExpansionStage>());
        p.run(ctx);

        if (ctx.args.landscape_only) {
            std::cout << "Landscape-only run complete; no selected primer panel was produced.\n";
            return 0;
        }
        {
        RunLogBlock final_reports("Final panel reports");
        print_solution(ctx.solution_primers);

        annotate_solution_dimers(ctx.segments, ctx.solution_primers,
            ctx.args.dimer_dg_cutoff.value_or(Args::default_dimer_dg_cutoff),
            ctx.args.dimer_temp_min.value_or(ctx.args.temp),
            ctx.args.dimer_temp_max.value_or(ctx.args.temp), ctx.args.output_file);
        save_solution(ctx.segments, ctx.solution_primers, ctx.args.output_file + ".bed");
        save_coverage_report(ctx.segments, ctx.solution_primers, ctx.args.output_file);
        std::size_t affected = 0;
        for (const auto& pair : ctx.solution_primers)
            affected += !pair.left.dimer_risk.partners.empty() + !pair.right.dimer_risk.partners.empty();
        final_reports.complete(std::to_string(ctx.solution_primers.size()) + " amplicons; " +
            std::to_string(affected) + " affected primers; " + ctx.args.output_file + ".bed");
        }

        return 0;
    } catch (const std::exception& e) {
        run_log.archive_error_if_needed();
        std::cerr << "Error: " << e.what() << "\n\n";
        print_usage(argv[0]);
        run_log.finish(false, e.what());
        return 2;
    }
}
