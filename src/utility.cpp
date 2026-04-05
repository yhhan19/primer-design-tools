#include "utility.hpp"

std::size_t read_fasta(const std::string &file_name, 
                       std::vector<std::string> &labels, 
                       std::vector<std::string> &data) {
    std::ifstream fin(file_name);
    if (! fin.is_open()) {
        std::cout << "unable to open file " << file_name << std::endl;
        return 1;
    }
    labels.clear(); data.clear();
    std::string line;
    while (getline(fin, line)) {
        if (line[0] == '>') {
            labels.push_back(line.substr(1));
            data.push_back("");
        }
        else {
            for (std::size_t i = 0; i < line.size(); i ++) {
                if (line[i] >= 'a' && line[i] <= 'z')
                    line[i] = line[i] - 32;
            }
            data[data.size() - 1] += line;
        }
    }
    fin.close();

    std::size_t min, max = -1;
    for (std::size_t i = 0; i < labels.size(); i ++) {
        if (min == -1 || data[i].size() < min) 
            min = data[i].size();
        if (max == -1 || data[i].size() > max) 
            max = data[i].size();
    }
    
    return labels.size();
}

bool is_ws(char c) {
    return c == '\n' || c == '\r' || c == '\t' || c == ' ';
}

char comp(char c) {
    switch (c) {
        case 'A': return 'T';
        case 'T': return 'A';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'a': return 't';
        case 't': return 'a';
        case 'c': return 'g';
        case 'g': return 'c';
        default:  return 'N';
    }
}

bool is_gc(char c) {
    return c == 'G' || c == 'C';
}

std::string revcomp(const std::string& s) {
    std::string rc;
    for (auto it = s.rbegin(); it != s.rend(); ++it)
        rc.push_back(comp(*it));
    return rc;
}

std::unordered_map<std::string, std::vector<Result>>
group_by_label(const std::vector<Result>& results) {
    std::unordered_map<std::string, std::vector<Result>> grouped;
    grouped.reserve(results.size());

    // label -> set(ref_index already seen)
    std::unordered_map<std::string, std::unordered_set<std::size_t>> seen;
    seen.reserve(results.size());

    for (const auto& r : results) {
        auto& s = seen[r.label];
        if (s.insert(r.ref_index).second) {
            grouped[r.label].push_back(r);
        }
    }

    return grouped;
}

void write_results(const std::string& output_file,
                   const std::vector<Result>& results,
                   const std::vector<std::string>& labels) {
    std::ofstream out(output_file);
    if (!out)
        throw std::runtime_error("Cannot open file: " + output_file);

    auto grouped = group_by_label(results);
    
    std::size_t count = 0;

    for (const auto& label : labels) {
        auto it = grouped.find(label);
        if (it == grouped.end()) continue;

        auto& vec = it->second;

        std::sort(vec.begin(), vec.end(),
                [](const Result& a, const Result& b) {
                    if (a.flag != b.flag) return a.flag < b.flag;
                    if (a.dg   != b.dg)   return a.dg   < b.dg;
                    return a.ref_index < b.ref_index;
                });

        out << ">" << label << "\n";

        for (const auto& r : vec) {
            out << r.data << '\t'
                << r.ref_seq << '\t'
                << r.ref_label << '\t'
                << r.ref_index << '\t'
                << r.flag << '\t'
                << r.dg << '\n';
            count ++;
        }
        out << "\n";
    }

    std::cout << "Result size      : " << results.size() << "\n"
              << "Output size      : " << count          << "\n"
              << "Output saved in  : " << output_file    << "\n";
}


std::vector<risk_t> random_risk(std::size_t size) {
    std::vector<risk_t> input;
    for (std::size_t i = 0; i < size; i ++) 
        input.push_back(4 * (float)(rand()) / (float)(RAND_MAX));
    return input;
}

std::size_t random_between(std::size_t min, std::size_t max) {
    std::size_t r = rand() * RAND_MAX + rand();
    return min + r % (max - min);
}

key_t to_key(index_t f, index_t r, index_t f_, index_t r_) {
    return (key_t) f + ((key_t) r << 16) 
        + ((key_t) f_ << 32) + ((key_t) r_ << 48);
}

key_t to_key_2(index_t r, index_t r_) {
    return (key_t) r + ((key_t) (r - r_) << 16);
}

void to_index(key_t k, index_t &f, index_t &r, index_t &f_, index_t &r_) {
    f = k & 0xffff;
    r = (k >> 16) & 0xffff;
    f_ = (k >> 32) & 0xffff;
    r_ = (k >> 48) & 0xffff;
}

void to_index_2(key_t k, index_t &r, index_t &r_) {
    r = k & 0xffff;
    r_ = r - (k >> 16) & 0xffff;
}

std::size_t read_rate(std::string filename, std::vector<risk_t> &input) {
    std::ifstream fin(filename);
    if (! fin.is_open()) {
        std::cout << "file not open" << std::endl;
        return -1;
    }
    input.clear();
    std::string line;
    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        risk_t rate;
        ss >> rate >> rate;
        input.push_back(rate);
    }
    fin.close();
    return 0;
}

std::size_t read_rate_2(std::string filename, 
                        std::vector<risk_t> &input, 
                        std::string &ref) {
    std::ifstream fin(filename);
    if (! fin.is_open()) {
        std::cout << "file not open" << std::endl;
        return -1;
    }
    input.clear();
    std::string line;
    std::vector<std::pair<index_t, risk_t>> temp;
    index_t max = 0;
    std::getline(fin, line);
    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        index_t i;
        risk_t rate;
        char comma, nc;
        ss >> i >> comma 
           >> nc >> comma 
           >> rate >> comma 
           >> rate >> comma 
           >> rate >> comma 
           >> rate >> comma >> rate;
        temp.push_back(std::make_pair(i, rate));
        ref += nc;
        max = std::max(max, i);
    }
    for (index_t i = 0; i < max; i ++) 
        input.push_back(0);
    for (auto elem : temp) {
        input[elem.first] = elem.second;
    }
    fin.close();
    return 0;
}

std::size_t read_rate_ref(std::string filename, 
                          std::vector<risk_t> &input, 
                          std::string &ref) {
    std::ifstream fin(filename);
    if (! fin.is_open()) {
        std::cout << "file not open" << std::endl;
        return 1;
    }
    input.clear();
    std::string line;
    std::vector<std::pair<index_t, risk_t>> temp;
    index_t max = 0;
    std::getline(fin, line);
    while (std::getline(fin, line)) {
        std::stringstream ss(line);
        index_t i;
        risk_t rate;
        char comma, nc;
        ss >> i >> comma >> nc >> comma >> rate;
        temp.push_back(std::make_pair(i, rate));
        ref += nc;
        max = std::max(max, i);
    }
    for (index_t i = 0; i < max; i ++) 
        input.push_back(0);
    for (auto elem : temp) {
        input[elem.first] = elem.second;
    }
    fin.close();
    return 0;
}

std::string display_PDR(std::vector<index_t> PDR) {
    std::string PDRs = "";
    for (index_t i = 0; i < PDR.size(); i ++) 
        PDRs += std::to_string(PDR[i]) + " ";
    return PDRs;
}

std::string msa_consensus(const std::vector<std::string>& msa) {
    char gap = '-';
    if (msa.empty()) return {};

    const std::size_t width = msa[0].size();
    for (const auto& seq : msa)
        if (seq.size() != width)
            throw std::invalid_argument(
                "msa_consensus: sequences differ in length"
            );

    std::string consensus;
    consensus.reserve(width);

    for (std::size_t col = 0; col < width; ++col) {
        std::array<int, ALPHABET_SIZE> freq{};
        int gap_count = 0;
        int total     = 0;

        for (const auto& seq : msa) {
            unsigned char c = static_cast<unsigned char>(seq[col]);
            if (seq[col] == gap) {
                ++gap_count;
            } else {
                ++freq[c];
                ++total;
            }
        }

        if (total == 0) {
            consensus += gap;
            continue;
        }

        char best   = gap;
        int  best_f = gap_count;

        for (int i = 0; i < 128; ++i) {
            if (best == gap || freq[i] > best_f) {
                best_f = freq[i];
                best   = static_cast<char>(i);
            }
        }

        consensus += best;
    }

    return consensus;
}

void display_primer_output(const PrimerOutput& out) {

    std::cout << out.pairs.size() << "\t"
              << out.left_oligos.size() << "\t"
              << out.right_oligos.size() << "\n";

    return ;
    auto print_oligo = [](const Oligo& o, const std::string& label) {
        std::cout << "  " << label
                  << "  seq="   << o.seq
                  << "  start=" << o.start
                  << "  len="   << o.length
                  << "  Tm="    << std::fixed << std::setprecision(2) << o.tm
                  << "  GC="    << std::fixed << std::setprecision(2) << o.gc
                  << "\n";
    };

    std::cout << "=== Pairs (" << out.pairs.size() << ") ===\n";
    for (std::size_t i = 0; i < out.pairs.size(); i++) {
        const PrimerResult& p = out.pairs[i];
        std::cout << "Pair " << i << "  product=" << p.product_size << "\n";
        print_oligo(p.left,  "L");
        print_oligo(p.right, "R");
    }

    std::cout << "\n=== Left oligos (" << out.left_oligos.size() << ") ===\n";
    for (std::size_t i = 0; i < out.left_oligos.size(); i++) {
        std::cout << "[" << i << "]";
        print_oligo(out.left_oligos[i], "");
    }

    std::cout << "\n=== Right oligos (" << out.right_oligos.size() << ") ===\n";
    for (std::size_t i = 0; i < out.right_oligos.size(); i++) {
        std::cout << "[" << i << "]";
        print_oligo(out.right_oligos[i], "");
    }
}

void convert(const std::vector<PrimerOutput>& primer_outputs,
                       std::vector<std::string>& labels,
                       std::vector<std::string>& sequences) {

    labels.clear();
    sequences.clear();
    
    for (size_t i = 0; i < primer_outputs.size(); ++i) {
        const auto& output = primer_outputs[i];
        
        // Extract primer pairs
        for (size_t j = 0; j < output.pairs.size(); ++j) {
            // Forward primer
            labels.push_back("OUT" + std::to_string(i) + "_P" + std::to_string(j) + "_F");
            sequences.push_back(output.pairs[j].left.seq);
            
            // Reverse primer
            labels.push_back("OUT" + std::to_string(i) + "_P" + std::to_string(j) + "_R");
            sequences.push_back(output.pairs[j].right.seq);
        }

        // Extract standalone left oligos
        for (size_t j = 0; j < output.left_oligos.size(); ++j) {
            labels.push_back("OUT" + std::to_string(i) + "_L" + std::to_string(j));
            sequences.push_back(output.left_oligos[j].seq);
        }
        
        // Extract standalone right oligos
        for (size_t j = 0; j < output.right_oligos.size(); ++j) {
            labels.push_back("OUT" + std::to_string(i) + "_R" + std::to_string(j));
            sequences.push_back(output.right_oligos[j].seq);
        }
    }
}

Oligo extract_oligo(const primer_rec& p,
                    const char* sequence,
                    bool is_right) {
    Oligo o;
    o.length = p.length;
    o.tm     = p.temp;
    o.gc     = p.gc_content;
    o.start  = is_right ? p.start - p.length + 1 : p.start;
    o.seq    = std::string(sequence + o.start, p.length);
    return o;
}

PrimerOutput extract_all(const p3retval* retval,
                         const seq_args*  sa) {
    PrimerOutput out;

    out.pairs.reserve(retval->best_pairs.num_pairs);
    for (int i = 0; i < retval->best_pairs.num_pairs; i++) {
        const primer_pair& pp = retval->best_pairs.pairs[i];
        PrimerResult r;
        r.left         = extract_oligo(*pp.left,  sa->sequence, false);
        r.right        = extract_oligo(*pp.right, sa->sequence, true);
        r.product_size = pp.product_size;
        out.pairs.push_back(r);
    }

    out.left_oligos.reserve(retval->fwd.num_elem);
    for (int i = 0; i < retval->fwd.num_elem; i++)
        out.left_oligos.push_back(
            extract_oligo(retval->fwd.oligo[i], sa->sequence, false));

    out.right_oligos.reserve(retval->rev.num_elem);
    for (int i = 0; i < retval->rev.num_elem; i++)
        out.right_oligos.push_back(
            extract_oligo(retval->rev.oligo[i], sa->sequence, true));

    return out;
}

std::vector<PrimerOutput> filterByDG(const std::vector<Result>& results,
                                     double dg_threshold,
                                     const std::vector<PrimerOutput>& original_primers) {
    
    // Group results by output index - these are primers to REMOVE
    std::unordered_map<int, std::vector<const Result*>> results_by_output;
    
    std::cout << "=== DG-BASED PRIMER REMOVAL ===" << std::endl;
    std::cout << "dG threshold: " << dg_threshold << std::endl;
    
    size_t total_to_remove = 0;
    for (const auto& result : results) {
        if (result.dg <= dg_threshold) {  // These meet criteria for REMOVAL
            total_to_remove++;
            if (result.label.length() >= 3 && result.label.substr(0, 3) == "OUT") {
                size_t pos = result.label.find('_');
                if (pos != std::string::npos) {
                    try {
                        int output_idx = std::stoi(result.label.substr(3, pos - 3));
                        if (output_idx >= 0 && output_idx < static_cast<int>(original_primers.size())) {
                            results_by_output[output_idx].push_back(&result);
                        }
                    } catch (...) {
                        continue;
                    }
                }
            }
        }
    }
    
    std::cout << "Total primers to remove: " << total_to_remove << "/" << results.size() 
              << " (dG <= " << dg_threshold << ")" << std::endl;
    
    std::vector<PrimerOutput> filtered_outputs(original_primers.size());
    
    // Collect and deduplicate removed primers
    if (total_to_remove > 0) {
        // Track the worst (lowest) dG for each unique primer
        std::map<std::tuple<int, std::string, int>, const Result*> worst_pairs;      // output, type(P), index -> result
        std::map<std::tuple<int, std::string, int>, const Result*> worst_standalone; // output, type(L/R), index -> result
        
        for (const auto& [output_idx, output_results] : results_by_output) {
            for (const auto* result : output_results) {
                const std::string& label = result->label;
                
                // Parse and categorize
                if (label.find("_P") != std::string::npos) {
                    // Pair primers
                    size_t p_pos = label.find("_P");
                    size_t f_pos = label.find("_F");
                    size_t r_pos = label.find("_R");
                    
                    if (p_pos != std::string::npos && (f_pos != std::string::npos || r_pos != std::string::npos)) {
                        try {
                            std::string idx_str = label.substr(p_pos + 2, 
                                (f_pos != std::string::npos ? f_pos : r_pos) - p_pos - 2);
                            int pair_idx = std::stoi(idx_str);
                            std::string direction = (f_pos != std::string::npos) ? "F" : "R";
                            
                            auto key = std::make_tuple(output_idx, "P" + std::to_string(pair_idx) + "_" + direction, pair_idx);
                            
                            // Keep the one with lowest (most negative) dG
                            if (!worst_pairs.count(key) || result->dg < worst_pairs[key]->dg) {
                                worst_pairs[key] = result;
                            }
                        } catch (...) {}
                    }
                }
                else if (label.find("_L") != std::string::npos || label.find("_R") != std::string::npos) {
                    // Standalone primers
                    char type = (label.find("_L") != std::string::npos) ? 'L' : 'R';
                    size_t type_pos = label.find(std::string("_") + type);
                    
                    try {
                        std::string idx_str = label.substr(type_pos + 2);
                        int idx = std::stoi(idx_str);
                        
                        auto key = std::make_tuple(output_idx, std::string(1, type), idx);
                        
                        // Keep the one with lowest (most negative) dG
                        if (!worst_standalone.count(key) || result->dg < worst_standalone[key]->dg) {
                            worst_standalone[key] = result;
                        }
                    } catch (...) {}
                }
            }
        }
        
        // TABLE A1: Filtered Primer Pairs
        if (!worst_pairs.empty()) {
            std::cout << "\nFiltered Primer Pairs" << std::endl;
            std::cout << std::string(73, '-') << std::endl;
            std::cout << std::setw(8) << "PDR pair"
                      << std::setw(8) << "Index"
                      << std::setw(10) << "Direction"
                      << std::setw(12) << "dG"
                      << std::setw(35) << "Sequence"
                      << std::endl;
            std::cout << std::string(73, '-') << std::endl;
            
            // Sort by output, then index, then direction
            std::vector<std::pair<std::tuple<int, std::string, int>, const Result*>> sorted_pairs(
                worst_pairs.begin(), worst_pairs.end());
            
            std::sort(sorted_pairs.begin(), sorted_pairs.end(), 
                [](const auto& a, const auto& b) {
                    if (std::get<0>(a.first) != std::get<0>(b.first)) 
                        return std::get<0>(a.first) < std::get<0>(b.first);
                    if (std::get<2>(a.first) != std::get<2>(b.first))
                        return std::get<2>(a.first) < std::get<2>(b.first);
                    return std::get<1>(a.first) < std::get<1>(b.first);
                });
            
            for (const auto& [key, result] : sorted_pairs) {
                int output_idx = std::get<0>(key);
                int pair_idx = std::get<2>(key);
                std::string type_dir = std::get<1>(key);
                
                std::string direction = type_dir.substr(type_dir.length() - 1);
                std::string reason = (direction == "F") ? "Forward primer" : "Reverse primer";
                
                std::string display_seq = result->data;
                if (display_seq.length() > 30) {
                    display_seq = display_seq.substr(0, 27) + "...";
                }
                
                std::cout << std::setw(8) << output_idx
                          << std::setw(8) << pair_idx
                          << std::setw(10) << direction
                          << std::setw(12) << std::fixed << std::setprecision(1) << result->dg
                          << std::setw(35) << display_seq
                          << std::endl;
            }
            std::cout << std::string(73, '-') << std::endl;
        }
        
        // TABLE A2: Filtered Standalone Primers
        if (!worst_standalone.empty()) {
            std::cout << "\nFiltered Oligos" << std::endl;
            std::cout << std::string(73, '-') << std::endl;
            std::cout << std::setw(8) << "PDR pair"
                      << std::setw(8) << "Index"
                      << std::setw(10) << "Type"
                      << std::setw(12) << "dG"
                      << std::setw(35) << "Sequence"
                      << std::endl;
            std::cout << std::string(73, '-') << std::endl;
            
            // Sort by output, then type, then index
            std::vector<std::pair<std::tuple<int, std::string, int>, const Result*>> sorted_standalone(
                worst_standalone.begin(), worst_standalone.end());
            
            std::sort(sorted_standalone.begin(), sorted_standalone.end(),
                [](const auto& a, const auto& b) {
                    if (std::get<0>(a.first) != std::get<0>(b.first))
                        return std::get<0>(a.first) < std::get<0>(b.first);
                    if (std::get<1>(a.first) != std::get<1>(b.first))
                        return std::get<1>(a.first) < std::get<1>(b.first);
                    return std::get<2>(a.first) < std::get<2>(b.first);
                });
            
            for (const auto& [key, result] : sorted_standalone) {
                int output_idx = std::get<0>(key);
                std::string type = std::get<1>(key);
                int idx = std::get<2>(key);
                
                std::string reason = (type == "L") ? "Standalone left" : "Standalone right";
                
                std::string display_seq = result->data;
                if (display_seq.length() > 30) {
                    display_seq = display_seq.substr(0, 27) + "...";
                }
                
                std::cout << std::setw(8) << output_idx
                          << std::setw(8) << idx
                          << std::setw(10) << type
                          << std::setw(12) << std::fixed << std::setprecision(1) << result->dg
                          << std::setw(35) << display_seq
                          << std::endl;
            }
            std::cout << std::string(73, '-') << std::endl;
        }
    }
    
    // Process each output for summary (same logic as before but simplified)
    struct OutputSummary {
        size_t original_pairs, original_left, original_right;
        size_t final_pairs, final_left, final_right;
        size_t removed_count;
        bool has_removal;
    };
    
    std::vector<OutputSummary> output_summaries;
    size_t total_pairs_kept = 0, total_left_kept = 0, total_right_kept = 0;
    
    for (size_t output_idx = 0; output_idx < original_primers.size(); ++output_idx) {
        const auto& original = original_primers[output_idx];
        PrimerOutput filtered_output;
        
        OutputSummary summary;
        summary.original_pairs = original.pairs.size();
        summary.original_left = original.left_oligos.size();
        summary.original_right = original.right_oligos.size();
        summary.removed_count = results_by_output.count(output_idx) ? results_by_output[output_idx].size() : 0;
        summary.has_removal = (summary.removed_count > 0);
        
        // Apply same filtering logic as before...
        std::set<int> remove_pair_indices_F, remove_pair_indices_R;
        std::set<int> remove_left_indices, remove_right_indices;
        
        if (results_by_output.count(output_idx)) {
            for (const auto* result : results_by_output[output_idx]) {
                const std::string& label = result->label;
                
                if (label.find("_P") != std::string::npos) {
                    size_t p_pos = label.find("_P");
                    size_t f_pos = label.find("_F");
                    size_t r_pos = label.find("_R");
                    
                    if (p_pos != std::string::npos && (f_pos != std::string::npos || r_pos != std::string::npos)) {
                        try {
                            std::string idx_str = label.substr(p_pos + 2, 
                                (f_pos != std::string::npos ? f_pos : r_pos) - p_pos - 2);
                            int pair_idx = std::stoi(idx_str);
                            
                            if (f_pos != std::string::npos) {
                                remove_pair_indices_F.insert(pair_idx);
                            } else {
                                remove_pair_indices_R.insert(pair_idx);
                            }
                        } catch (...) {}
                    }
                }
                else if (label.find("_L") != std::string::npos) {
                    size_t l_pos = label.find("_L");
                    try {
                        int left_idx = std::stoi(label.substr(l_pos + 2));
                        remove_left_indices.insert(left_idx);
                    } catch (...) {}
                }
                else if (label.find("_R") != std::string::npos) {
                    size_t r_pos = label.find("_R");
                    try {
                        int right_idx = std::stoi(label.substr(r_pos + 2));
                        remove_right_indices.insert(right_idx);
                    } catch (...) {}
                }
            }
        }
        
        // Apply filtering and count
        summary.final_pairs = 0;
        for (size_t i = 0; i < original.pairs.size(); ++i) {
            int pair_idx = static_cast<int>(i);
            if (!remove_pair_indices_F.count(pair_idx) && !remove_pair_indices_R.count(pair_idx)) {
                filtered_output.pairs.push_back(original.pairs[i]);
                summary.final_pairs++;
            }
        }
        
        summary.final_left = 0;
        for (size_t i = 0; i < original.left_oligos.size(); ++i) {
            int left_idx = static_cast<int>(i);
            if (!remove_left_indices.count(left_idx)) {
                filtered_output.left_oligos.push_back(original.left_oligos[i]);
                summary.final_left++;
            }
        }
        
        summary.final_right = 0;
        for (size_t i = 0; i < original.right_oligos.size(); ++i) {
            int right_idx = static_cast<int>(i);
            if (!remove_right_indices.count(right_idx)) {
                filtered_output.right_oligos.push_back(original.right_oligos[i]);
                summary.final_right++;
            }
        }
        
        total_pairs_kept += summary.final_pairs;
        total_left_kept += summary.final_left;
        total_right_kept += summary.final_right;
        
        output_summaries.push_back(summary);
        filtered_outputs[output_idx] = std::move(filtered_output);
    }
    
    // TABLE B: Summary by output
    std::cout << "\nOutput Summary" << std::endl;
    std::cout << std::string(62, '-') << std::endl;
    std::cout << std::setw(8) << "PDR pair"
              << std::setw(18) << "Original (P/L/R)"
              << std::setw(18) << "Removed (P/L/R)"
              << std::setw(18) << "Final (P/L/R)"
              << std::endl;
    std::cout << std::string(62, '-') << std::endl;
    
    for (size_t i = 0; i < output_summaries.size(); ++i) {
        const auto& summary = output_summaries[i];
        
        std::string original_str = std::to_string(summary.original_pairs) + "/" + 
                                  std::to_string(summary.original_left) + "/" + 
                                  std::to_string(summary.original_right);
        
        std::string final_str = std::to_string(summary.final_pairs) + "/" + 
                               std::to_string(summary.final_left) + "/" + 
                               std::to_string(summary.final_right);
        
        // Calculate actual removed counts per type
        size_t removed_pairs = summary.original_pairs - summary.final_pairs;
        size_t removed_left = summary.original_left - summary.final_left;
        size_t removed_right = summary.original_right - summary.final_right;
        
        std::string removed_str = std::to_string(removed_pairs) + "/" + 
                                 std::to_string(removed_left) + "/" + 
                                 std::to_string(removed_right);
        
        std::string status = summary.has_removal ? "Filtered" : "No removal";
        
        std::cout << std::setw(8) << i
                  << std::setw(18) << original_str
                  << std::setw(18) << removed_str
                  << std::setw(18) << final_str
                  << std::endl;
    }
    
    // Calculate totals
    size_t total_original = 0, total_final = 0;
    for (const auto& summary : output_summaries) {
        total_original += summary.original_pairs + summary.original_left + summary.original_right;
    }
    total_final = total_pairs_kept + total_left_kept + total_right_kept;
    
    std::cout << std::string(62, '-') << std::endl;
    std::cout << std::setw(8) << "TOTAL"
              << std::setw(18) << std::to_string(total_original)
              << std::setw(18) << total_to_remove
              << std::setw(18) << std::to_string(total_final)
              << std::endl;
    std::cout << std::string(62, '-') << std::endl;
    
    return filtered_outputs;
}
