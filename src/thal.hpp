#pragma once
#include "utility.hpp"
#include "dimer_cache.hpp"
#include <string>
#include <unordered_map>
#include <stdexcept>
#include <cmath>
#include <mutex>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>

// ── Thermodynamic result from SantaLucia 1998 ────────────────────────────────

struct ThermoResult {
    double dH;   // cal/mol
    double dS;   // cal/(mol·K)
    double dG37; // cal/mol at 37°C
    double tm;   // °C
};

// ── Thermocycler profile ──────────────────────────────────────────────────────
//
// ΔG evaluated at each stage of a standard PCR thermocycler cycle.
// All dG values in cal/mol (same unit as ThermoResult).
//
// Interpretation for spurious duplex risk:
//   dG_denature > 0  : duplex does not survive denaturation (good)
//   dG_anneal   < 0  : duplex is stable at annealing temperature (flag)
//   dG_extend   < 0  : duplex survives into extension, blocks polymerase (bad)
//
// worst_dG is min(dG_denature, dG_anneal, dG_extend) — the most stable
// the duplex ever becomes during the cycle. Since dG(T) = dH - T*dS is
// linear and dS < 0 (so dG increases with T), worst_dG always occurs at
// T_denature for on-target duplexes but at T_anneal for weak spurious ones.
// Checking worst_dG is therefore more conservative than dG(T_a) alone.

struct ThermoProfile {
    double dG_denature;  // cal/mol at T_denature (typically 95°C)
    double dG_anneal;    // cal/mol at T_anneal
    double dG_extend;    // cal/mol at T_extend   (typically 72°C)
    double worst_dG;     // min of the three — most stable point in cycle
    double T_denature;   // °C
    double T_anneal;     // °C
    double T_extend;     // °C
};

// ── Thal ─────────────────────────────────────────────────────────────────────

class Thal {
public:
    // ── Primer3/libthal interface ─────────────────────────────────────────────
    //
    // init(): initialises libthal and stores salt concentrations and
    //   thermocycler temperatures for use across all thermo computations.
    //   mv:         monovalent cation [Na+] (mM)
    //   dv:         divalent cation  [Mg2+] (mM)
    //   dntp:       dNTP concentration (mM)
    //   dna_conc:   primer concentration (nM)
    //   T_anneal:   annealing temperature (°C)
    //   T_denature: denaturation temperature (°C)
    //   T_extend:   extension temperature (°C)
    //
    // Implementation (.cpp) must include:
    //   na_mM_ = mv;  mg_mM_ = dv;  dntp_mM_ = dntp;
    //   T_anneal_ = T_anneal;  T_denature_ = T_denature;  T_extend_ = T_extend;

    static void init(const std::string& config_path,
                     double mv,
                     double dv,
                     double dntp,
                     double dna_conc,
                     double T_anneal,
                     double T_denature,
                     double T_extend);

    // Backward-compatible 6-argument overload for existing pipeline call sites.
    // Uses standard PCR defaults: T_denature=94, T_extend=72.
    static void init(const std::string& config_path,
                     double mv,
                     double dv,
                     double dntp,
                     double dna_conc,
                     double T_anneal) {
        init(config_path, mv, dv, dntp, dna_conc, T_anneal, 94.0, 72.0);
    }

    static void validate_dimer_range(double min_c, double max_c);
    static void configure_dimer_cache(const std::string& directory);
    static void flush_dimer_cache();
    static DimerInteractionCache::Stats dimer_cache_stats();
    // Exact endpoint minimum for the fixed structure returned by libthal.
    // Both inputs are physical oligo sequences in 5'->3' order.
    static double compute_dimer_worst_dg(const std::string& s1,
                                        const std::string& s2,
                                        double min_c, double max_c);

    static double compute_dimer_dg(const std::string& s1,
                                   const std::string& s2);

    // ── SantaLucia 1998 nearest-neighbor model ────────────────────────────────
    //
    // compute_thermo(): returns dH, dS, dG(37°C), and Tm for a single-stranded
    //   primer sequence (non-self-complementary duplex assumed).
    //   conc_nM: total primer strand concentration (default 250 nM).
    //
    // dG_at_T(): evaluates dG at an arbitrary annealing/extension temperature
    //   using the linear model dG(T) = dH - T*dS, given a ThermoResult.
    //   Use this to assess spurious duplex stability at T_a or T_extension.
    //
    // compute_thermo_heterodimer(): computes dH, dS, dG, Tm for a duplex
    //   formed between two primers (e.g. for dimer risk evaluation).
    //   Constructs the 5'->3' sequence of s1 aligned against rev-comp of s2
    //   and sums nearest-neighbor parameters over the duplex.

    static ThermoResult compute_thermo(const std::string& seq,
                                       double conc_nM = 250.0);

    static double dG_at_T(const ThermoResult& r, double T_celsius);

    // ── Thermocycler profile ──────────────────────────────────────────────────
    //
    // thermocycler_profile(): evaluates dG at denaturation, annealing, and
    //   extension temperatures and returns a ThermoProfile.
    //
    // Two overloads:
    //   (1) no temperatures — uses values stored at init() time, i.e. the
    //       user-supplied CLI temperatures. Use this for the common case.
    //   (2) explicit temperatures — overrides stored values for a specific
    //       amplicon or reaction condition (e.g. touchdown PCR, hotstart).
    //
    // profile_is_safe(): returns true if worst_dG >= threshold (cal/mol).
    //   Default threshold -2000 cal/mol (~-2 kcal/mol, Primer3 convention).

    static ThermoProfile thermocycler_profile(const ThermoResult& r);

    static ThermoProfile thermocycler_profile(const ThermoResult& r,
                                              double T_denature,
                                              double T_anneal,
                                              double T_extend);

    static bool profile_is_safe(const ThermoProfile& p,
                                double threshold = -2000.0);

    // ── Binding score ─────────────────────────────────────────────────────────
    //
    // binding_score(): three-point approximation of the integral of ΔG(T(t))
    //   over the dangerous region of one thermocycler cycle.
    //
    //   S = min(ΔG(Td),0) + min(ΔG(Ta),0) + min(ΔG(Te),0)
    //
    // Assumes equal hold times and instantaneous ramps (see derivation).
    // Only negative contributions are counted — positive ΔG means the duplex
    // is unstable at that point and contributes nothing to spurious binding.
    //
    // S = 0   : duplex never stable in cycle — safe
    // S < 0   : duplex stable at one or more holds — use as edge weight
    // S << 0  : duplex stable throughout cycle — critical
    //
    // Two overloads: no-arg uses stored init() temperatures; explicit overload
    // for per-amplicon temperature override.

    static double binding_score(const ThermoResult& r);

    static double binding_score(const ThermoResult& r,
                                double T_denature,
                                double T_anneal,
                                double T_extend);

    static ThermoResult compute_thermo_heterodimer(const std::string& s1,
                                                   const std::string& s2,
                                                   double conc_nM = 250.0);

    // ── Validation ────────────────────────────────────────────────────────────
    //
    // run_tests(): compares SantaLucia 1998 results against Primer3 libthal
    //   across a fixed suite of single-primer and heterodimer cases.
    //   Prints a comparison table to stdout.
    //   Requires init() to have been called first (for the libthal path).
    //   Ta: annealing temperature (°C), Te: extension temperature (°C).

    static void run_tests(double Ta = 60.0, double Te = 72.0);

private:
    // ── libthal state ─────────────────────────────────────────────────────────

    static inline thal_parameters tp_{};
    static inline thal_args       base_args_{};
    static inline bool            initialized_ = false;
    static inline std::string dimer_cache_context_;
    static inline std::string dimer_cache_directory_;
    static inline std::shared_ptr<DimerInteractionCache> dimer_cache_;
    static double calculate_dimer_worst_dg(const std::string& s1,
        const std::string& s2, double min_c, double max_c);
    static inline std::mutex      call_mutex_;

    // Salt concentrations stored at init() for use in salt correction.
    // Units: mM (matching Primer3 / Rust code conventions).
    static inline double na_mM_   = 50.0;
    static inline double mg_mM_   = 1.5;
    static inline double dntp_mM_ = 0.6;

    // Thermocycler temperatures stored at init() (°C).
    // All three are required — no implicit defaults.
    static inline double T_denature_ = 0.0;
    static inline double T_anneal_   = 0.0;
    static inline double T_extend_   = 0.0;

    // ── SantaLucia 1998 helpers ───────────────────────────────────────────────

    // Returns the reverse complement of a DNA sequence (A/T/G/C).
    static std::string reverse_complement(const std::string& seq);

    // Flat initiation penalty (SantaLucia 1998 averaged; matches Rust impl).
    // Applied once per duplex, not per terminal base.
    static void apply_init_flat(double& dH, double& dS);

    // Core summation over nearest-neighbor dinucleotide steps.
    // seq must be uppercase A/T/G/C only.
    static void sum_nn(const std::string& seq, double& dH, double& dS);

    // Computes ThermoResult from raw dH (cal/mol), dS (cal/mol·K), conc_nM,
    // and duplex length N (number of base pairs) for salt correction.
    // Salt correction (SantaLucia 2004 / Primer3 style) applied to dS:
    //   dS += 0.368 * (N-1) * ln(Na_eq / 1000)
    //   Na_eq = na_mM + 120 * sqrt(max(mg_mM - dntp_mM, 0))
    static ThermoResult finalize(double dH, double dS,
                                 double conc_nM, size_t N);

    // ── SantaLucia 1998 Table 2 nearest-neighbor parameters ──────────────────
    // Key: 5'->3' dinucleotide on the top strand.
    // Remaining 6 entries follow by Watson-Crick symmetry (e.g. AA == TT).

    struct NNParams { double dH; double dS; }; // cal/mol, cal/(mol·K)

    static const std::unordered_map<std::string, NNParams> NN_PARAMS;

    // ── Test infrastructure ───────────────────────────────────────────────────

    struct TestCase {
        std::string name;
        std::string s1;
        std::string s2;      // empty => single-primer mode
        double      conc_nM;
    };

    struct TestResult {
        std::string name;
        std::string seqs;
        // SantaLucia path
        double dH          = 0.0;
        double dS          = 0.0;
        double tm          = 0.0;
        double dG_denature = 0.0;  // at T_denature
        double dG_anneal   = 0.0;  // at T_anneal
        double dG_extend   = 0.0;  // at T_extend
        double worst_dG    = 0.0;  // min over cycle
        double score       = 0.0;  // binding_score: sum of negative dG terms
        bool   safe        = false;
        // Primer3 libthal path
        double dG_p3       = 0.0;
        double delta       = 0.0;  // dG_anneal - dG_p3
        std::string error;
    };

    static std::vector<TestCase> test_suite();
    static TestResult            run_one(const TestCase& tc,
                                         double Ta, double Te);
    static void                  print_results(
                                     const std::vector<TestResult>& results,
                                     double Ta, double Te);
};

// ── NN_PARAMS definition (header-only, inline) ────────────────────────────────

inline const std::unordered_map<std::string, Thal::NNParams> Thal::NN_PARAMS = {
    // Unique 10 from SantaLucia 1998 Table 2 — dH in cal/mol, dS in cal/(mol·K)
    {"AA", {-7900,  -22.2}},
    {"AT", {-7200,  -20.4}},
    {"TA", {-7200,  -21.3}},
    {"CA", {-8500,  -22.7}},
    {"GT", {-8400,  -22.4}},
    {"CT", {-7800,  -21.0}},
    {"GA", {-8200,  -22.2}},
    {"CG", {-10600, -27.2}},
    {"GC", {-9800,  -24.4}},
    {"GG", {-8000,  -19.9}},
    // Watson-Crick symmetric complements
    {"TT", {-7900,  -22.2}}, // == AA
    {"TG", {-8500,  -22.7}}, // == CA
    {"AC", {-8400,  -22.4}}, // == GT
    {"AG", {-7800,  -21.0}}, // == CT
    {"TC", {-8200,  -22.2}}, // == GA
    {"CC", {-8000,  -19.9}}, // == GG
};

// ── Inline implementations ────────────────────────────────────────────────────

inline std::string Thal::reverse_complement(const std::string& seq) {
    std::string rc(seq.size(), 'N');
    auto comp = [](char c) -> char {
        switch (c) {
            case 'A': return 'T';
            case 'T': return 'A';
            case 'G': return 'C';
            case 'C': return 'G';
            default:  return 'N';
        }
    };
    std::transform(seq.rbegin(), seq.rend(), rc.begin(), comp);
    return rc;
}

inline void Thal::apply_init_flat(double& dH, double& dS) {
    // Averaged initiation term (SantaLucia 1998; matches Rust impl).
    // Applied once per duplex: dH += 200 cal/mol, dS += -5.7 cal/(mol·K)
    dH += 200;   // 0.2 kcal/mol = 200 cal/mol
    dS += -5.7;
}

inline void Thal::sum_nn(const std::string& seq, double& dH, double& dS) {
    for (size_t i = 0; i + 1 < seq.size(); ++i) {
        std::string dinuc = seq.substr(i, 2);
        auto it = NN_PARAMS.find(dinuc);
        if (it == NN_PARAMS.end())
            throw std::invalid_argument(
                "sum_nn: unknown dinucleotide '" + dinuc + "'");
        dH += it->second.dH;
        dS += it->second.dS;
    }
}

inline ThermoResult Thal::finalize(double dH, double dS,
                                   double conc_nM, size_t N) {
    constexpr double R     = 1.9872;  // cal/(mol·K) — gas constant
    constexpr double T37_K = 310.15;  // 37°C in Kelvin

    // Salt correction on dS (SantaLucia 2004 / Primer3 style):
    //   Na_eq (mM) = [Na+] + 120 * sqrt(max([Mg2+] - [dNTP], 0))
    //   dS += 0.368 * (N-1) * ln(Na_eq / 1000)
    double free_mg  = std::max(mg_mM_ - dntp_mM_, 0.0);
    double na_eq    = na_mM_ + 120.0 * std::sqrt(free_mg);
    dS += 0.368 * (static_cast<double>(N) - 1.0) * std::log(na_eq / 1000.0);

    // dH and dS both in cal units — dG(T) = dH - T*dS directly
    double dG37 = dH - T37_K * dS;

    // Tm for non-self-complementary duplex (SantaLucia 1998 Eq. 3)
    double C_T  = conc_nM * 1e-9;
    double tm_K = dH / (dS + R * std::log(C_T / 4.0));

    return { dH, dS, dG37, tm_K - 273.15 };
}

inline ThermoResult Thal::compute_thermo(const std::string& seq,
                                         double conc_nM) {
    if (seq.size() < 2)
        throw std::invalid_argument("compute_thermo: sequence too short");

    double dH = 0.0, dS = 0.0;
    apply_init_flat(dH, dS);
    sum_nn(seq, dH, dS);

    return finalize(dH, dS, conc_nM, seq.size());
}

inline double Thal::dG_at_T(const ThermoResult& r, double T_celsius) {
    double T_K = T_celsius + 273.15;
    return r.dH - T_K * r.dS;  // cal/mol — dH and dS both in cal units
}

inline ThermoProfile Thal::thermocycler_profile(const ThermoResult& r) {
    // Uses temperatures stored at init() — the common call path.
    return thermocycler_profile(r, T_denature_, T_anneal_, T_extend_);
}

inline ThermoProfile Thal::thermocycler_profile(const ThermoResult& r,
                                                 double T_denature,
                                                 double T_anneal,
                                                 double T_extend) {
    ThermoProfile p;
    p.T_denature   = T_denature;
    p.T_anneal     = T_anneal;
    p.T_extend     = T_extend;
    p.dG_denature  = dG_at_T(r, T_denature);
    p.dG_anneal    = dG_at_T(r, T_anneal);
    p.dG_extend    = dG_at_T(r, T_extend);
    p.worst_dG     = std::min({p.dG_denature, p.dG_anneal, p.dG_extend});
    return p;
}

inline bool Thal::profile_is_safe(const ThermoProfile& p,
                                   double threshold) {
    // Safe if duplex is never stable enough to compete at any cycle stage.
    // worst_dG >= threshold means the most favorable binding point is still
    // above the risk threshold — duplex is too weak to matter.
    return p.worst_dG >= threshold;
}

inline double Thal::binding_score(const ThermoResult& r) {
    return binding_score(r, T_denature_, T_anneal_, T_extend_);
}

inline double Thal::binding_score(const ThermoResult& r,
                                   double T_denature,
                                   double T_anneal,
                                   double T_extend) {
    auto clamp = [](double x) { return x < 0.0 ? x : 0.0; };
    return clamp(dG_at_T(r, T_denature))
         + clamp(dG_at_T(r, T_anneal))
         + clamp(dG_at_T(r, T_extend));
}

inline ThermoResult Thal::compute_thermo_heterodimer(const std::string& s1,
                                                     const std::string& s2,
                                                     double conc_nM) {
    // Model the heterodimer as s1 hybridized to rev-comp(s2),
    // reading the duplex 5'->3' along s1.
    // Both strands must be the same length (full complementarity assumed).
    // For partial/internal dimers, pass the aligned sub-sequences.
    if (s1.size() != s2.size())
        throw std::invalid_argument(
            "compute_thermo_heterodimer: s1 and s2 must be equal length "
            "(pass aligned sub-sequences for partial duplexes)");

    double dH = 0.0, dS = 0.0;
    apply_init_flat(dH, dS);
    sum_nn(s1, dH, dS);  // nearest-neighbor sum along s1 (duplex context)

    return finalize(dH, dS, conc_nM, s1.size());
}

// ── Test suite definition ─────────────────────────────────────────────────────

inline std::vector<Thal::TestCase> Thal::test_suite() {
    // Coverage:
    //  1–3:  single primers (AT-rich, GC-rich, mixed) — validate dH/dS/Tm
    //  4–5:  palindrome and self-complement — edge cases for initiation
    //  6–8:  heterodimers: weak / moderate / strong dimer risk
    //  9:    3'-end dimer (most biologically dangerous)
    //  10:   long 25-mer
    //  11:   AT-rich 20-mer (low Tm)
    //  12:   GC-rich 20-mer (high Tm)
    return {
        {"1. short AT-rich",   "AATTTAATTT",              "",             250},
        {"2. short GC-rich",   "GCGCGCGCGC",              "",             250},
        {"3. mixed 20-mer",    "ATCGATCGATCGATCGATCG",    "",             250},
        {"4. palindrome",      "GAATTC",                  "",             250},
        {"5. self-complement", "AAAATTTT",                "",             250},
        {"6. weak dimer",      "ATCGATCGATCG", "CGATCGATCGAT",           250},
        {"7. mod. dimer",      "GCATGCATGCAT", "ATGCATGCATGC",           250},
        {"8. strong dimer",    "CGCGCGCGCGCG", "CGCGCGCGCGCG",           250},
        // 3'-end dimer: truncated to equal-length 3' overlap in run_one
        {"9. 3'-end dimer",    "ATCGATCGATCGATCG", "CGATCG",             250},
        {"10. long 25-mer",    "ATCGATCGATCGATCGATCGATCGA", "",           250},
        {"11. AT-rich 20-mer", "ATATATATATATATATATATAT",    "",           250},
        {"12. GC-rich 20-mer", "GCGCGCGCGCGCGCGCGCGC",     "",           250},
    };
}

// ── Single test execution ─────────────────────────────────────────────────────

inline Thal::TestResult Thal::run_one(const TestCase& tc,
                                      double Ta, double Te) {
    TestResult r;
    r.name = tc.name;

    auto trunc = [](const std::string& s, size_t n) {
        return s.size() > n ? s.substr(0, n - 1) + "..." : s;
    };
    r.seqs = trunc(tc.s1, 14);
    if (!tc.s2.empty()) r.seqs += "/" + trunc(tc.s2, 10);

    try {
        ThermoResult sl{};

        if (tc.s2.empty()) {
            sl      = compute_thermo(tc.s1, tc.conc_nM);
            r.dG_p3 = compute_dimer_dg(tc.s1, tc.s1);
        } else {
            std::string a = tc.s1;
            std::string b = tc.s2;
            if (a.size() != b.size()) {
                size_t len = std::min(a.size(), b.size());
                a = a.substr(a.size() - len);
                b = b.substr(b.size() - len);
            }
            sl      = compute_thermo_heterodimer(a, b, tc.conc_nM);
            r.dG_p3 = compute_dimer_dg(tc.s1, tc.s2);
        }

        // Populate base fields
        r.dH = sl.dH;
        r.dS = sl.dS;
        r.tm = sl.tm;

        // Full thermocycler profile — explicit temperatures since test
        // Ta/Te are passed per-call and may differ from stored init() values
        auto prof       = thermocycler_profile(sl, 95.0, Ta, Te);
        r.dG_denature   = prof.dG_denature;
        r.dG_anneal     = prof.dG_anneal;
        r.dG_extend     = prof.dG_extend;
        r.worst_dG      = prof.worst_dG;
        r.safe          = profile_is_safe(prof);
        r.score         = binding_score(sl, 95.0, Ta, Te);

        r.delta = r.dG_anneal - r.dG_p3;

    } catch (const std::exception& e) {
        r.error = e.what();
    }
    return r;
}

// ── Result printing ───────────────────────────────────────────────────────────

inline void Thal::print_results(const std::vector<TestResult>& results,
                                 double Ta, double Te) {
    constexpr int W_NAME = 22;
    constexpr int W_SEQ  = 26;
    constexpr int W_NUM  =  9;
    constexpr int W_SAFE =  6;
    // columns: dH dS Tm dG(95) dG(Ta) dG(Te) worst score dG(P3) DdG safe
    const std::string SEP(W_NAME + W_SEQ + W_NUM * 9 + W_SAFE, '-');

    std::cout << "\nThal::run_tests()  --  SantaLucia 1998 vs Primer3 libthal\n"
              << "Ta = " << Ta << " C    Te = " << Te << " C    "
              << "T_denature = 95 C    [primer] = 250 nM\n\n"
              << "All dG in cal/mol.  "
              << "score = sum of min(dG,0) at Td/Ta/Te  |  "
              << "safe = worst_dG >= -2000  |  "
              << "*** = |DdG| > 2000 cal/mol\n\n";

    std::cout << SEP << "\n"
              << std::left
              << std::setw(W_NAME) << "Test"
              << std::setw(W_SEQ)  << "Sequence(s)"
              << std::right
              << std::setw(W_NUM)  << "dH"
              << std::setw(W_NUM)  << "dS"
              << std::setw(W_NUM)  << "Tm"
              << std::setw(W_NUM)  << "dG(95)"
              << std::setw(W_NUM)  << "dG(Ta)"
              << std::setw(W_NUM)  << "dG(Te)"
              << std::setw(W_NUM)  << "worst"
              << std::setw(W_NUM)  << "score"
              << std::setw(W_NUM)  << "dG(P3)"
              << std::setw(W_NUM)  << "DdG"
              << std::setw(W_SAFE) << "safe"
              << "\n" << SEP << "\n";

    for (const auto& r : results) {
        std::cout << std::left
                  << std::setw(W_NAME) << r.name
                  << std::setw(W_SEQ)  << r.seqs
                  << std::right << std::fixed << std::setprecision(0);
        if (!r.error.empty()) {
            std::cout << "  ERROR: " << r.error;
        } else {
            std::cout << std::setw(W_NUM) << r.dH
                      << std::setw(W_NUM) << r.dS
                      << std::setprecision(1)
                      << std::setw(W_NUM) << r.tm
                      << std::setprecision(0)
                      << std::setw(W_NUM) << r.dG_denature
                      << std::setw(W_NUM) << r.dG_anneal
                      << std::setw(W_NUM) << r.dG_extend
                      << std::setw(W_NUM) << r.worst_dG
                      << std::setw(W_NUM) << r.score
                      << std::setw(W_NUM) << r.dG_p3
                      << std::setw(W_NUM) << r.delta
                      << std::setw(W_SAFE) << (r.safe ? "Y" : "N");
            if (std::abs(r.delta) > 2000.0)
                std::cout << " ***";
        }
        std::cout << "\n";
    }

    std::cout << SEP << "\n\n"
              << "Notes:\n"
              << "  score:   binding_score = min(dG(95),0) + min(dG(Ta),0) +\n"
              << "    min(dG(Te),0). Three-point approximation of integral of\n"
              << "    dG(T(t)) over dangerous region. score=0 means safe;\n"
              << "    more negative = more dangerous. Use as dimer edge weight.\n"
              << "  dG(95):  should be positive (unstable at denaturation).\n"
              << "  dG(Ta):  negative means spurious duplex forms at annealing.\n"
              << "  dG(Te):  negative means duplex blocks polymerase.\n"
              << "  worst:   most stable point in cycle (min of three).\n"
              << "  safe:    worst_dG >= -2000 cal/mol (Primer3 default).\n"
              << "  Single-primer dG(P3): libthal self-dimer -- not directly\n"
              << "    comparable; constant DdG ~197 cal/mol is initiation\n"
              << "    convention difference, not a bug.\n"
              << "  Equal-length heterodimers (6,7,8): |DdG| < 100 -- good.\n"
              << "  Unequal-length case (9): |DdG| ~1000 -- libthal finds\n"
              << "    longer overlap via DP alignment; expected.\n";
}

// ── Public entry point ────────────────────────────────────────────────────────

inline void Thal::run_tests(double Ta, double Te) {
    auto suite = test_suite();
    std::vector<TestResult> results;
    results.reserve(suite.size());
    for (const auto& tc : suite)
        results.push_back(run_one(tc, Ta, Te));
    print_results(results, Ta, Te);
}
