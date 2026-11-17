#pragma once
#include "utility.hpp"
#include <functional>

// Energies are keyed by ordered physical 5'->3' sequences and temperature bounds.
// Graph indices, pool assignment, geometry, objectives and cutoffs are not keys.
class DimerInteractionCache {
public:
    struct Stats { std::size_t reused = 0, computed = 0; };
    DimerInteractionCache(std::string context, const std::string& directory);
    double get(const std::string& first, const std::string& second,
               double min_c, double max_c, const std::function<double()>& calculate);
    void flush();
    Stats stats() const;
private:
    std::string context_, directory_;
    mutable std::mutex mutex_;
    std::mutex flush_mutex_;
    std::unordered_map<std::string, std::shared_future<double>> entries_;
    std::map<std::string, double> pending_;
    Stats stats_;
};
