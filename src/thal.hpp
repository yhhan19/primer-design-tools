#pragma once
#include "utility.hpp"

class Thal {
public:
    static void init(const std::string& config_path,
                     double mv,
                     double dv,
                     double dntp,
                     double dna_conc,
                     double temp);
    static double compute_dimer_dg(const std::string& s1,
                                   const std::string& s2);

private:
    static inline thal_parameters tp_{};
    static inline thal_args base_args_{};
    static inline bool initialized_ = false;
    static inline std::mutex call_mutex_;
};
