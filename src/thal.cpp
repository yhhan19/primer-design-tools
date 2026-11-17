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

    // Bump the scoring revision when libthal or the range-score algorithm
    // changes. Store actual loaded parameter contents, not just their path.
    std::ostringstream context;
    context << "dpro-libthal-range-v1"
#ifdef DPRO_DIMER_ENGINE_ID
            << '|' << DPRO_DIMER_ENGINE_ID
#else
            << "|custom-build|" << __DATE__ << '|' << __TIME__
#endif
            << std::setprecision(std::numeric_limits<double>::max_digits10)
            << '|' << base_args_.type << '|' << base_args_.maxLoop
            << '|' << mv << '|' << dv << '|' << dntp << '|' << dna_conc
            << '|' << base_args_.dimer;
    for (const char* parameter : {tp_.dangle_dh, tp_.dangle_ds, tp_.loops_dh, tp_.loops_ds,
            tp_.stack_dh, tp_.stack_ds, tp_.stackmm_dh, tp_.stackmm_ds,
            tp_.tetraloop_dh, tp_.tetraloop_ds, tp_.triloop_dh, tp_.triloop_ds,
            tp_.tstack_tm_inf_ds, tp_.tstack_dh, tp_.tstack2_dh, tp_.tstack2_ds}) {
        context << '|';
        if (!parameter) { context << "null"; continue; }
        for (const unsigned char* c = reinterpret_cast<const unsigned char*>(parameter); *c; ++c)
            context << std::hex << std::setw(2) << std::setfill('0') << unsigned(*c);
    }
    dimer_cache_context_ = context.str();

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
    if (!dimer_cache_) return calculate_dimer_worst_dg(s1, s2, min_c, max_c);
    return dimer_cache_->get(s1, s2, min_c, max_c, [&] {
        return calculate_dimer_worst_dg(s1, s2, min_c, max_c);
    });
}

void Thal::configure_dimer_cache(const std::string& directory) {
    if (!initialized_) throw std::runtime_error("Thal not initialized. Call Thal::init()");
    if (directory == "off") {
        flush_dimer_cache();
        dimer_cache_.reset();
        dimer_cache_directory_.clear();
        return;
    }
    const auto path = fs::absolute(directory).lexically_normal().string();
    if (dimer_cache_ && path == dimer_cache_directory_) return;
    flush_dimer_cache();
    auto cache = std::make_shared<DimerInteractionCache>(dimer_cache_context_, path);
    dimer_cache_ = std::move(cache);
    dimer_cache_directory_ = path;
}

void Thal::flush_dimer_cache() {
    if (dimer_cache_) dimer_cache_->flush();
}

DimerInteractionCache::Stats Thal::dimer_cache_stats() {
    return dimer_cache_ ? dimer_cache_->stats() : DimerInteractionCache::Stats{};
}

double Thal::calculate_dimer_worst_dg(const std::string& s1,
    const std::string& s2, double min_c, double max_c) {
    validate_dimer_range(min_c, max_c);
    if (!initialized_) throw std::runtime_error("Thal not initialized. Call Thal::init()");
    auto args = base_args_; // Never mutate shared thermodynamic arguments.
    args.temp = min_c - KELVIN_ZERO;
    thal_results r{};
    thal((const unsigned char*)s1.c_str(),
         (const unsigned char*)s2.c_str(), &args, THL_FAST, &r);
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
