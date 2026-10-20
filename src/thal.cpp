#include "thal.hpp"

void Thal::init(const std::string& config_path,
                double mv,
                double dv,
                double dntp,
                double dna_conc,
                double T_anneal,
                double T_denature,
                double T_extend) {

    if (initialized_) return;

    if (thal_set_null_parameters(&tp_) != 0) {
        throw std::runtime_error("thal_set_null_parameters failed");
    }

    thal_results o{};
    if (thal_load_parameters(config_path.c_str(), &tp_, &o) != 0) {
        throw std::runtime_error(
            std::string("thal_load_parameters failed: ") + o.msg
        );
    }

    set_thal_oligo_default_args(&base_args_);
    base_args_.type     = thal_any;
    base_args_.mv       = mv;
    base_args_.dv       = dv;
    base_args_.dntp     = dntp;
    base_args_.dna_conc = dna_conc;
    base_args_.temp     = T_anneal - (KELVIN_ZERO);  // libthal uses Kelvin
    base_args_.dimer    = 1;

    // Store salt and temperatures for SantaLucia path
    na_mM_      = mv;
    mg_mM_      = dv;
    dntp_mM_    = dntp;
    T_anneal_   = T_anneal;
    T_denature_ = T_denature;
    T_extend_   = T_extend;

    initialized_ = true;
}

double Thal::compute_dimer_dg(const std::string& s1,
                              const std::string& s2) {
    if (!initialized_) {
        throw std::runtime_error("Thal not initialized. Call Thal::init()");
    }

    thal_results r{};

    // std::lock_guard<std::mutex> lk(call_mutex_);

    thal((const unsigned char*)s1.c_str(),
         (const unsigned char*)revcomp(s2).c_str(),
         &base_args_,
         THL_FAST,
         &r);

    if (r.temp == THAL_ERROR_SCORE) {
        throw std::runtime_error(std::string("thal error: ") + r.msg);
    }

    return r.dg;
}
void Thal::validate_dimer_range(double min_c, double max_c) {
    if (!std::isfinite(min_c) || !std::isfinite(max_c) ||
        min_c <= KELVIN_ZERO || max_c <= KELVIN_ZERO || min_c > max_c) {
        throw std::invalid_argument(
            "Dimer temperatures must be finite, above absolute zero, and min <= max");
    }
}

double Thal::compute_dimer_worst_dg(const std::string& s1,
                                   const std::string& s2,
                                   double min_c, double max_c) {
    validate_dimer_range(min_c, max_c);
    if (!initialized_) throw std::runtime_error("Thal not initialized. Call Thal::init()");
    auto args = base_args_; // Never mutate shared thermodynamic arguments.
    args.temp = min_c - KELVIN_ZERO;
    thal_results r{};
    const auto s2_rc = revcomp(s2); // Preserve existing sequence convention.
    thal((const unsigned char*)s1.c_str(),
         (const unsigned char*)s2_rc.c_str(), &args, THL_FAST, &r);
    if (r.temp == THAL_ERROR_SCORE) {
        throw std::runtime_error(std::string("thal error: ") + r.msg);
    }
    if (!std::isfinite(r.dh) || !std::isfinite(r.ds)) {
        throw std::runtime_error("Non-finite dimer thermodynamics");
    }
    // libthal returns salt-corrected dS in cal/(mol K), dH in cal/mol.
    // Its dimer structure search uses a fixed reference temperature; this is
    // a fixed-structure range score, not an ensemble or a structure re-search.
    // No structure leaves the zero-initialized dH/dS at zero (legacy behavior).
    const double dg_min = r.dh - (min_c - KELVIN_ZERO) * r.ds;
    const double dg_max = r.dh - (max_c - KELVIN_ZERO) * r.ds;
    return std::min(dg_min, dg_max);
}
