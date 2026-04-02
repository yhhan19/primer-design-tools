#include "risk_optimizer.hpp"

double RiskOptimizer::calculate_nucleotide_diversity(const std::vector<char>& bases) {
    if (bases.size() <= 1) return 0.0;
    
    std::unordered_map<char, int> counts;
    for (char base : bases) {
        counts[base]++;
    }
    
    int total_pairs = bases.size() * (bases.size() - 1) / 2;
    int different_pairs = 0;
    
    for (auto it1 = counts.begin(); it1 != counts.end(); ++it1) {
        for (auto it2 = std::next(it1); it2 != counts.end(); ++it2) {
            different_pairs += it1->second * it2->second;
        }
    }
    
    return static_cast<double>(different_pairs) / total_pairs;
}

double RiskOptimizer::calculate_gap_penalty(const std::vector<char>& position_data) {
    int gap_count = std::count(position_data.begin(), position_data.end(), '-');
    double gap_fraction = static_cast<double>(gap_count) / position_data.size();
    return std::sqrt(gap_fraction);  // Square root penalty
}

double RiskOptimizer::calculate_ambiguity_penalty(const std::vector<char>& position_data) {
    std::unordered_set<char> ambiguous_codes = {'N', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V'};
    
    int ambiguous_count = 0;
    for (char base : position_data) {
        if (ambiguous_codes.count(base) > 0) {
            ambiguous_count++;
        }
    }
    
    return static_cast<double>(ambiguous_count) / position_data.size();
}

double RiskOptimizer::calculate_shannon_entropy(const std::vector<char>& bases) {
    if (bases.empty()) return std::log2(4.0);  // Maximum entropy
    
    std::unordered_map<char, int> counts;
    for (char base : bases) {
        counts[base]++;
    }
    
    double entropy = 0.0;
    for (const auto& pair : counts) {
        if (pair.second > 0) {
            double p = static_cast<double>(pair.second) / bases.size();
            entropy -= p * std::log2(p);
        }
    }
    
    // Normalize by maximum possible entropy
    double max_entropy = std::log2(std::min(4.0, static_cast<double>(counts.size())));
    return max_entropy > 0 ? entropy / max_entropy : 0.0;
}

double RiskOptimizer::calculate_conservation_score(const std::vector<char>& bases) {
    if (bases.empty()) return 1.0;
    
    std::unordered_map<char, int> counts;
    for (char base : bases) {
        counts[base]++;
    }
    
    int max_count = 0;
    for (const auto& pair : counts) {
        max_count = std::max(max_count, pair.second);
    }
    
    double most_common_freq = static_cast<double>(max_count) / bases.size();
    return 1.0 - most_common_freq;
}

double RiskOptimizer::calculate_gc_content_deviation(const std::vector<std::string>& sequences, 
                                    size_t position, int window_size) {
    std::vector<double> gc_contents;
    
    for (const auto& seq : sequences) {
        // Define window boundaries
        int start = std::max(0, static_cast<int>(position) - window_size);
        int end = std::min(static_cast<int>(seq.length()), static_cast<int>(position) + window_size + 1);
        
        // Extract window and count valid bases
        std::vector<char> valid_bases;
        for (int pos = start; pos < end; ++pos) {
            char base = std::toupper(seq[pos]);
            if (base == 'A' || base == 'T' || base == 'G' || base == 'C') {
                valid_bases.push_back(base);
            }
        }
        
        if (valid_bases.empty()) {
            gc_contents.push_back(0.5);  // Neutral
            continue;
        }
        
        // Calculate GC content for this sequence's window
        int gc_count = 0;
        for (char base : valid_bases) {
            if (base == 'G' || base == 'C') {
                gc_count++;
            }
        }
        
        double gc_content = static_cast<double>(gc_count) / valid_bases.size();
        gc_contents.push_back(gc_content);
    }
    
    if (gc_contents.empty()) return 1.0;
    
    // Average GC content across all sequences
    double avg_gc = std::accumulate(gc_contents.begin(), gc_contents.end(), 0.0) / gc_contents.size();
    
    // Calculate deviation from optimal 50%
    double deviation = std::abs(avg_gc - 0.5) / 0.5;
    return std::min(deviation, 1.0);
}

double RiskOptimizer::calculate_local_complexity(const std::vector<std::string>& sequences, 
                                size_t position, int window_size) {
    std::vector<double> complexity_scores;
    
    for (const auto& seq : sequences) {
        // Define window boundaries
        int start = std::max(0, static_cast<int>(position) - window_size / 2);
        int end = std::min(static_cast<int>(seq.length()), static_cast<int>(position) + window_size / 2 + 1);
        
        // Extract window without gaps
        std::string window;
        for (int pos = start; pos < end; ++pos) {
            char base = std::toupper(seq[pos]);
            if (base != '-') {
                window += base;
            }
        }
        
        if (window.length() < 3) {
            complexity_scores.push_back(1.0);  // Maximum risk for short windows
            continue;
        }
        
        // Count unique 2-mers (dinucleotides)
        std::unordered_set<std::string> unique_kmers;
        for (size_t i = 0; i < window.length() - 1; ++i) {
            unique_kmers.insert(window.substr(i, 2));
        }
        
        int max_possible_kmers = std::min(static_cast<int>(window.length() - 1), 16);  // 4^2 = 16
        double complexity = static_cast<double>(unique_kmers.size()) / max_possible_kmers;
        
        complexity_scores.push_back(1.0 - complexity);  // Convert to risk
    }
    
    return complexity_scores.empty() ? 1.0 : 
            std::accumulate(complexity_scores.begin(), complexity_scores.end(), 0.0) / complexity_scores.size();
}

double RiskOptimizer::calculate_indel_risk(const std::vector<std::string>& sequences, 
                            size_t position, int window_size) {
    int indel_count = 0;
    
    for (const auto& seq : sequences) {
        // Define window boundaries
        int start = std::max(0, static_cast<int>(position) - window_size);
        int end = std::min(static_cast<int>(seq.length()), static_cast<int>(position) + window_size + 1);
        
        // Count gaps in window
        for (int pos = start; pos < end; ++pos) {
            if (seq[pos] == '-') {
                indel_count++;
            }
        }
    }
    
    int total_positions = sequences.size() * (2 * window_size + 1);
    double indel_frequency = static_cast<double>(indel_count) / total_positions;
    
    return std::min(indel_frequency * 2.0, 1.0);  // Scale and cap at 1.0
}

double RiskOptimizer::calculate_tm_stability_risk(const std::vector<char>& valid_bases, 
                                    const std::vector<std::string>& sequences,
                                    size_t position, int window_size) {
    // Simplified Tm calculation based on nearest neighbor approximation
    // Higher variance in Tm across sequences = higher risk
    
    std::vector<double> tm_estimates;
    
    for (const auto& seq : sequences) {
        int start = std::max(0, static_cast<int>(position) - window_size / 2);
        int end = std::min(static_cast<int>(seq.length()), static_cast<int>(position) + window_size / 2 + 1);
        
        std::string window;
        for (int pos = start; pos < end; ++pos) {
            char base = std::toupper(seq[pos]);
            if (base != '-' && base != 'N') {
                window += base;
            }
        }
        
        if (window.length() < 10) {
            tm_estimates.push_back(50.0);  // Default Tm
            continue;
        }
        
        // Simplified Tm calculation (Wallace rule + GC correction)
        int at_count = 0, gc_count = 0;
        for (char base : window) {
            if (base == 'A' || base == 'T') at_count++;
            else if (base == 'G' || base == 'C') gc_count++;
        }
        
        double tm = 64.9 + 41.0 * (gc_count - 16.4) / window.length();
        if (window.length() < 14) {
            tm = (at_count * 2) + (gc_count * 4);  // Wallace rule for short oligos
        }
        
        tm_estimates.push_back(tm);
    }
    
    if (tm_estimates.empty()) return 1.0;
    
    // Calculate variance in Tm estimates
    double mean_tm = std::accumulate(tm_estimates.begin(), tm_estimates.end(), 0.0) / tm_estimates.size();
    double variance = 0.0;
    
    for (double tm : tm_estimates) {
        variance += (tm - mean_tm) * (tm - mean_tm);
    }
    variance /= tm_estimates.size();
    
    double std_dev = std::sqrt(variance);
    
    // Risk increases with Tm variance (primers should have similar Tm)
    // Normalize: >5°C difference is high risk
    return std::min(std_dev / 5.0, 1.0);
}

double RiskOptimizer::calculate_secondary_structure_risk(const std::vector<char>& valid_bases,
                                        const std::vector<std::string>& sequences,
                                        size_t position, int window_size) {
    // Detect potential secondary structures (hairpins, self-complementarity)
    std::vector<double> structure_risks;
    
    for (const auto& seq : sequences) {
        int start = std::max(0, static_cast<int>(position) - window_size / 2);
        int end = std::min(static_cast<int>(seq.length()), static_cast<int>(position) + window_size / 2 + 1);
        
        std::string window;
        for (int pos = start; pos < end; ++pos) {
            char base = std::toupper(seq[pos]);
            if (base != '-' && base != 'N') {
                window += base;
            }
        }
        
        if (window.length() < 8) {
            structure_risks.push_back(0.0);
            continue;
        }
        
        double risk = 0.0;
        
        // Check for simple hairpin potential (palindromic sequences)
        int hairpin_matches = 0;
        int check_length = std::min(window.length() / 2, 6UL);
        
        for (size_t i = 0; i < check_length; ++i) {
            char forward = window[i];
            char reverse = window[window.length() - 1 - i];
            
            // Check for Watson-Crick complementarity
            bool complementary = false;
            if ((forward == 'A' && reverse == 'T') || (forward == 'T' && reverse == 'A') ||
                (forward == 'G' && reverse == 'C') || (forward == 'C' && reverse == 'G')) {
                complementary = true;
            }
            
            if (complementary) {
                hairpin_matches++;
            }
        }
        
        double hairpin_risk = static_cast<double>(hairpin_matches) / check_length;
        
        // Check for homopolymer runs (AAAA, GGGG, etc.)
        int max_run_length = 1;
        int current_run = 1;
        
        for (size_t i = 1; i < window.length(); ++i) {
            if (window[i] == window[i-1]) {
                current_run++;
                max_run_length = std::max(max_run_length, current_run);
            } else {
                current_run = 1;
            }
        }
        
        double homopolymer_risk = 0.0;
        if (max_run_length >= 4) {
            homopolymer_risk = std::min((max_run_length - 3) / 5.0, 1.0);
        }
        
        risk = std::max(hairpin_risk, homopolymer_risk);
        structure_risks.push_back(risk);
    }
    
    return structure_risks.empty() ? 0.0 :
            std::accumulate(structure_risks.begin(), structure_risks.end(), 0.0) / structure_risks.size();
}

 
std::vector<risk_t> RiskOptimizer::compute_risk(const std::vector<std::string>& sequences) {
    if (sequences.empty()) return {};
    
    size_t length = sequences[0].size();
    std::vector<risk_t> risk_scores(length, 0.0);
    
    // Configuration parameters
    const int gc_window_size = 10;     // Window for GC content analysis
    const int complexity_window_size = 10;  // Window for complexity analysis
    const int indel_window_size = 5;   // Window for indel risk analysis
    const int tm_window_size = 20;     // Window for Tm stability analysis
    const int structure_window_size = 15; // Window for secondary structure analysis
    
    RiskWeights weights;
    
    std::cout << "Computing comprehensive risk scores for " << length << " positions..." << std::endl;
    
    for (size_t pos = 0; pos < length; ++pos) {
        if (pos % 1000 == 0) {
            std::cout << "  Position " << pos << "/" << length << std::endl;
        }
        
        // Extract bases at this position
        std::vector<char> position_data;
        std::vector<char> valid_bases;
        
        for (const auto& seq : sequences) {
            char base = std::toupper(seq[pos]);
            position_data.push_back(base);
            if (base != '-' && base != 'N') {
                valid_bases.push_back(base);
            }
        }
        
        if (valid_bases.empty()) {
            risk_scores[pos] = 1.0;  // Maximum risk
            continue;
        }

        /*
        struct RiskWeights {
            double diversity = 0.25;
            double gaps = 0.20;
            double ambiguity = 0.15;
            double entropy = 0.15;
            double conservation = 0.10;
            double gc_deviation = 0.05;
            double complexity = 0.05;
            double indels = 0.05;
        };*/
        
        // Calculate all risk components
        double diversity_risk = calculate_nucleotide_diversity(valid_bases);
        double gap_risk = calculate_gap_penalty(position_data);
        double ambiguity_risk = calculate_ambiguity_penalty(position_data);
        double entropy_risk = calculate_shannon_entropy(valid_bases);
        double conservation_risk = calculate_conservation_score(valid_bases);
        double gc_risk = calculate_gc_content_deviation(sequences, pos, gc_window_size);
        double complexity_risk = calculate_local_complexity(sequences, pos, complexity_window_size);
        double indel_risk = calculate_indel_risk(sequences, pos, indel_window_size);
        
        // Additional advanced risk components
        double tm_risk = calculate_tm_stability_risk(valid_bases, sequences, pos, tm_window_size);
        double structure_risk = calculate_secondary_structure_risk(valid_bases, sequences, pos, structure_window_size);
        
        // Adjust weights to include new components
        double total_weight = weights.diversity + weights.gaps + weights.ambiguity + 
                            weights.entropy + weights.conservation + weights.gc_deviation +
                            weights.complexity + weights.indels + 0.05 + 0.05; // tm + structure
        
        // Weighted combination of all risk factors
        double total_risk = (
            weights.diversity * diversity_risk +
            weights.gaps * gap_risk +
            weights.ambiguity * ambiguity_risk +
            weights.entropy * entropy_risk +
            weights.conservation * conservation_risk +
            weights.gc_deviation * gc_risk +
            weights.complexity * complexity_risk +
            weights.indels * indel_risk +
            0.05 * tm_risk +           // Tm stability weight
            0.05 * structure_risk      // Secondary structure weight
        ) / total_weight;
        
        risk_scores[pos] = std::min(total_risk, 1.0);
    }
    
    std::cout << "Risk computation completed!" << std::endl;
    return risk_scores;
}