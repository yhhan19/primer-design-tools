#include "thal.hpp"

void Thal::init(const std::string& config_path,
                double mv,
                double dv,
                double dntp,
                double dna_conc,
                double temp) {

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
    base_args_.type = thal_any;
    base_args_.mv = mv;
    base_args_.dv = dv;
    base_args_.dntp = dntp;
    base_args_.dna_conc = dna_conc;
    base_args_.temp = temp - (KELVIN_ZERO);
    base_args_.dimer = 1;

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
