#include "utility.hpp"
#include "thal.hpp"

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
    static const std::string valid = "ATCGNRYSWKMBDHV-";
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
            for (char c : line) {
                if (valid.find(c) != std::string::npos) {
                    data[data.size() - 1] += c;
                }
            }
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

void write_fasta(const std::string& filename,
                 const std::vector<std::string>& headers,
                 const std::vector<std::string>& sequences,
                 int line_width) {
    if (headers.size() != sequences.size()) {
        throw std::invalid_argument("Headers and sequences must have the same length");
    }

    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    for (size_t i = 0; i < headers.size(); i++) {
        file << ">" << headers[i] << "\n";
        for (size_t j = 0; j < sequences[i].size(); j += line_width) {
            file << sequences[i].substr(j, line_width) << "\n";
        }
    }
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

    // std::cout << "Output saved in  : " << output_file    << "\n";
}


std::vector<risk_t> random_risk(std::size_t size) {
    std::vector<risk_t> input;
    for (std::size_t i = 0; i < size; i ++) 
        input.push_back(4 * (float)(rand()) / (float)(RAND_MAX));
    return input;
}

std::size_t random_between(std::size_t min, std::size_t max) {
    if (min == max) return min;
    std::size_t r = rand() * RAND_MAX + rand();
    return min + r % (max - min + 1);  // +1 to make range inclusive
}

key_t to_key(index_t f, index_t r, index_t f_, index_t r_) {
    return (key_t) f + ((key_t) r << 16) 
        + ((key_t) f_ << 32) + ((key_t) r_ << 48);
}

key_t to_key_2(index_t r, index_t r_) {
    return (static_cast<key_t>(r) << 32) | r_;
}

void to_index(key_t k, index_t &f, index_t &r, index_t &f_, index_t &r_) {
    f = k & 0xffff;
    r = (k >> 16) & 0xffff;
    f_ = (k >> 32) & 0xffff;
    r_ = (k >> 48) & 0xffff;
}

void to_index_2(key_t k, index_t &r, index_t &r_) {
    r = static_cast<index_t>(k >> 32);
    r_ = static_cast<index_t>(k);
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

void char_freq(const std::string& seq) {
    std::map<char, int> freq;
    for (char c : seq) {
        freq[c]++;
    }
    for (const auto& [c, count] : freq) {
        std::cout << "'" << c << "' (ASCII " << (int)c << "): " 
                  << count << " (" << 100.0 * count / seq.size() << "%)\n";
    }
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
        for (size_t j = 0; !output.pairs_required && j < output.pairs.size(); ++j) {
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

RemovalSets computeRemoval(const std::vector<const Result*>& output_results,
                           double threshold,
                           std::size_t num_return) {
    RemovalSets rs;
    rs.num_return = num_return;
    for (const auto* result : output_results) {
        if (result->dg > threshold) continue;
        const std::string& label = result->label;
        if (label.find("_P") != std::string::npos) {
            size_t p_pos = label.find("_P"), f_pos = label.find("_F"), r_pos = label.find("_R");
            if (p_pos != std::string::npos && (f_pos != std::string::npos || r_pos != std::string::npos)) {
                try {
                    int pair_idx = std::stoi(label.substr(p_pos + 2,
                        (f_pos != std::string::npos ? f_pos : r_pos) - p_pos - 2));
                    if (f_pos != std::string::npos) rs.remove_pair_F.insert(pair_idx);
                    else                             rs.remove_pair_R.insert(pair_idx);
                } catch (...) {}
            }
        } else if (label.find("_L") != std::string::npos) {
            try { rs.remove_left.insert(std::stoi(label.substr(label.find("_L") + 2))); } catch (...) {}
        } else if (label.find("_R") != std::string::npos) {
            try { rs.remove_right.insert(std::stoi(label.substr(label.find("_R") + 2))); } catch (...) {}
        }
    }
    return rs;
}

std::vector<PrimerOutput> filterByDG_relax(
    const std::vector<Result>& results,
    double dg_threshold,
    const std::vector<PrimerOutput>& original_primers,
    std::size_t num_return) {

    const double RELAX_STEP    = 500;
    const int    MAX_RELAX_STEPS = 10;

    // Group ALL results by output index (not just those meeting threshold)
    std::unordered_map<int, std::vector<const Result*>> results_by_output;
    for (const auto& result : results) {
        if (result.label.length() >= 3 && result.label.substr(0, 3) == "OUT") {
            size_t pos = result.label.find('_');
            if (pos != std::string::npos) {
                try {
                    int output_idx = std::stoi(result.label.substr(3, pos - 3));
                    if (output_idx >= 0 && output_idx < static_cast<int>(original_primers.size()))
                        results_by_output[output_idx].push_back(&result);
                } catch (...) {}
            }
        }
    }

    // Count total to remove at the global threshold (for header line)
    size_t total_to_remove = 0;
    for (const auto& result : results)
        if (result.dg <= dg_threshold) total_to_remove++;

    std::cout << "Total primers to remove: " << total_to_remove << "/" << results.size()
              << " (dG <= " << dg_threshold << ")" << std::endl;

    std::vector<PrimerOutput> filtered_outputs(original_primers.size());

    // Collect worst-dG per unique primer for reporting tables
    std::map<std::tuple<int, std::string, int>, const Result*> worst_pairs;
    std::map<std::tuple<int, std::string, int>, const Result*> worst_standalone;

    // Per-output effective thresholds (may differ from global if relaxed)
    std::vector<double> effective_threshold(original_primers.size(), dg_threshold);

    struct OutputSummary {
        size_t original_pairs, original_left, original_right;
        size_t final_pairs, final_left, final_right;
        bool has_removal;
        bool was_relaxed;
    };

    std::vector<OutputSummary> output_summaries;
    size_t total_pairs_kept = 0, total_left_kept = 0, total_right_kept = 0;

    for (size_t output_idx = 0; output_idx < original_primers.size(); ++output_idx) {
        const auto& original = original_primers[output_idx];
        PrimerOutput filtered_output;

        OutputSummary summary;
        summary.original_pairs = original.pairs.size();
        summary.original_left  = original.left_oligos.size();
        summary.original_right = original.right_oligos.size();
        summary.has_removal    = false;
        summary.was_relaxed    = false;

        auto& output_results = results_by_output[output_idx];  // empty vec if not in map

        double local_threshold = dg_threshold;
        RemovalSets rs = computeRemoval(output_results, local_threshold, num_return);

        if (rs.all_removed(original)) {
            for (int step = 1; step <= MAX_RELAX_STEPS; ++step) {
                local_threshold -= RELAX_STEP;
                rs = computeRemoval(output_results, local_threshold, num_return);
                if (!rs.all_removed(original)) {
                    std::cout << "[PDR " << output_idx << "] threshold relaxed "
                              << dg_threshold << " -> " << local_threshold << " kcal/mol\n";
                    effective_threshold[output_idx] = local_threshold;
                    summary.was_relaxed = true;
                    break;
                }
            }
            if (rs.all_removed(original))
                std::cerr << "[PDR " << output_idx << "] WARNING: still empty after max relaxation ("
                          << (dg_threshold - MAX_RELAX_STEPS * RELAX_STEP) << " kcal/mol)\n";
        }

        summary.has_removal = (!rs.remove_pair_F.empty() || !rs.remove_pair_R.empty() ||
                               !rs.remove_left.empty()   || !rs.remove_right.empty());

        // Collect worst-dG entries for reporting tables (at effective threshold)
        for (const auto* result : output_results) {
            if (result->dg > local_threshold) continue;
            const std::string& label = result->label;

            if (label.find("_P") != std::string::npos) {
                size_t p_pos = label.find("_P"), f_pos = label.find("_F"), r_pos = label.find("_R");
                if (p_pos != std::string::npos && (f_pos != std::string::npos || r_pos != std::string::npos)) {
                    try {
                        int pair_idx = std::stoi(label.substr(p_pos + 2,
                            (f_pos != std::string::npos ? f_pos : r_pos) - p_pos - 2));
                        std::string direction = (f_pos != std::string::npos) ? "F" : "R";
                        auto key = std::make_tuple((int)output_idx,
                            "P" + std::to_string(pair_idx) + "_" + direction, pair_idx);
                        if (!worst_pairs.count(key) || result->dg < worst_pairs[key]->dg)
                            worst_pairs[key] = result;
                    } catch (...) {}
                }
            } else if (label.find("_L") != std::string::npos || label.find("_R") != std::string::npos) {
                char type = (label.find("_L") != std::string::npos) ? 'L' : 'R';
                size_t type_pos = label.find(std::string("_") + type);
                try {
                    int idx = std::stoi(label.substr(type_pos + 2));
                    auto key = std::make_tuple((int)output_idx, std::string(1, type), idx);
                    if (!worst_standalone.count(key) || result->dg < worst_standalone[key]->dg)
                        worst_standalone[key] = result;
                } catch (...) {}
            }
        }

        // Apply filtering
        summary.final_pairs = 0;
        for (size_t i = 0; i < original.pairs.size(); ++i)
            if (!rs.remove_pair_F.count(i) && !rs.remove_pair_R.count(i)) {
                filtered_output.pairs.push_back(original.pairs[i]);
                summary.final_pairs++;
            }

        summary.final_left = 0;
        for (size_t i = 0; i < original.left_oligos.size(); ++i)
            if (!rs.remove_left.count(i)) {
                filtered_output.left_oligos.push_back(original.left_oligos[i]);
                summary.final_left++;
            }

        summary.final_right = 0;
        for (size_t i = 0; i < original.right_oligos.size(); ++i)
            if (!rs.remove_right.count(i)) {
                filtered_output.right_oligos.push_back(original.right_oligos[i]);
                summary.final_right++;
            }

        filtered_output.pairs_required = original.pairs_required;
        filtered_output.segment_id = original.segment_id;
        filtered_output.pdr_left = original.pdr_left;
        filtered_output.pdr_right = original.pdr_right;
        if (filtered_output.pairs_required) {
            retain_surviving_pairs(filtered_output);
            if (filtered_output.pairs.empty())
                throw std::runtime_error("Off-target filtering left no valid candidate pair for output " + std::to_string(output_idx));
            summary.final_pairs = filtered_output.pairs.size();
        }
        total_pairs_kept += summary.final_pairs;
        total_left_kept  += summary.final_left;
        total_right_kept += summary.final_right;

        output_summaries.push_back(summary);
        filtered_outputs[output_idx] = std::move(filtered_output);
    }

    // TABLE A1: Filtered Primer Pairs
    if (!worst_pairs.empty()) {
        /*
        std::cout << "\nFiltered Primer Pairs" << std::endl;
        std::cout << std::string(73, '-') << std::endl;
        std::cout << std::setw(8)  << "PDR pair"
                  << std::setw(8)  << "Index"
                  << std::setw(10) << "Direction"
                  << std::setw(12) << "dG"
                  << std::setw(35) << "Sequence"
                  << std::endl;
        std::cout << std::string(73, '-') << std::endl;
        */

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
            int output_idx   = std::get<0>(key);
            int pair_idx     = std::get<2>(key);
            std::string type_dir  = std::get<1>(key);
            std::string direction = type_dir.substr(type_dir.length() - 1);
            std::string display_seq = result->data;
            if (display_seq.length() > 30) display_seq = display_seq.substr(0, 27) + "...";
            /*
            std::cout << std::setw(8)  << output_idx
                      << std::setw(8)  << pair_idx
                      << std::setw(10) << direction
                      << std::setw(12) << std::fixed << std::setprecision(1) << result->dg
                      << std::setw(35) << display_seq
                      << std::endl;
            */
        }
        // std::cout << std::string(73, '-') << std::endl;
    }

    // TABLE A2: Filtered Standalone Primers
    if (!worst_standalone.empty()) {
        /*
        std::cout << "\nFiltered Oligos" << std::endl;
        std::cout << std::string(73, '-') << std::endl;
        std::cout << std::setw(8)  << "PDR pair"
                  << std::setw(8)  << "Index"
                  << std::setw(10) << "Type"
                  << std::setw(12) << "dG"
                  << std::setw(35) << "Sequence"
                  << std::endl;
        std::cout << std::string(73, '-') << std::endl;
        */

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
            int output_idx  = std::get<0>(key);
            std::string type = std::get<1>(key);
            int idx          = std::get<2>(key);
            std::string display_seq = result->data;
            if (display_seq.length() > 30) display_seq = display_seq.substr(0, 27) + "...";
            /*
            std::cout << std::setw(8)  << output_idx
                      << std::setw(8)  << idx
                      << std::setw(10) << type
                      << std::setw(12) << std::fixed << std::setprecision(1) << result->dg
                      << std::setw(35) << display_seq
                      << std::endl;
            */
        }
        // std::cout << std::string(73, '-') << std::endl;
    }

    // TABLE B: Summary
    std::cout << "\nSummary" << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    std::cout << std::setw(8)  << "PDR pair"
              << std::setw(18) << "Original (P/L/R)"
              << std::setw(18) << "Removed (P/L/R)"
              << std::setw(18) << "Final (P/L/R)"
              << std::setw(18) << "Threshold"
              << std::endl;
    std::cout << std::string(80, '-') << std::endl;

    size_t total_original = 0;
    for (size_t i = 0; i < output_summaries.size(); ++i) {
        const auto& summary = output_summaries[i];
        total_original += summary.original_pairs + summary.original_left + summary.original_right;

        std::string original_str = std::to_string(summary.original_pairs) + "/" +
                                   std::to_string(summary.original_left)  + "/" +
                                   std::to_string(summary.original_right);
        std::string final_str    = std::to_string(summary.final_pairs) + "/" +
                                   std::to_string(summary.final_left)  + "/" +
                                   std::to_string(summary.final_right);
        std::string removed_str  = std::to_string(summary.original_pairs - summary.final_pairs) + "/" +
                                   std::to_string(summary.original_left  - summary.final_left)  + "/" +
                                   std::to_string(summary.original_right - summary.final_right);
        std::string thresh_str   = std::to_string(effective_threshold[i]);
        if (summary.was_relaxed) thresh_str += "*"; else thresh_str += " ";

        std::cout << std::setw(8)  << i
                  << std::setw(18) << original_str
                  << std::setw(18) << removed_str
                  << std::setw(18) << final_str
                  << std::setw(18) << thresh_str
                  << std::endl;
    }

    size_t total_final = total_pairs_kept + total_left_kept + total_right_kept,
           total_remove = total_original - total_final;
    
    std::cout << std::string(80, '-') << std::endl;
    std::cout << std::setw(8)  << "TOTAL"
              << std::setw(18) << std::to_string(total_original)
              << std::setw(18) << std::to_string(total_remove)
              << std::setw(18) << std::to_string(total_final)
              << std::setw(18) << ""
              << std::endl;
    std::cout << std::string(80, '-') << std::endl;
    std::cout << "* threshold auto-relaxed to preserve at least one oligo per type" << std::endl;

    return filtered_outputs;
}

void print_solution(const std::vector<PrimerResult>& solution) {
    const int W_IDX  = 6;
    const int W_POOL = 6;
    const int W_SEG  = 6;
    const int W_PDR  = 16;
    const int W_SEQ  = 28;
    const int W_TM   = 7;
    const int W_GC   = 7;
    const int W_SIZE = 7;
    const int TOTAL  = W_IDX + W_POOL + W_SEG + W_PDR + W_SEQ*2 + W_TM*2 + W_GC*2 + W_SIZE + 2;

    std::cout << "\nSolution Primers\n"
              << std::string(TOTAL, '-') << "\n"
              << std::right
              << std::setw(W_IDX)  << "idx"
              << std::setw(W_POOL) << "pool"
              << std::setw(W_SEG)  << "seg"
              << std::setw(W_PDR)  << "PDR"
              << std::setw(W_SEQ)  << "left seq"
              << std::setw(W_TM)   << "Tm"
              << std::setw(W_GC)   << "GC%"
              << std::setw(W_SEQ)  << "right seq"
              << std::setw(W_TM)   << "Tm"
              << std::setw(W_GC)   << "GC%"
              << std::setw(W_SIZE) << "size"
              << "\n"
              << std::string(TOTAL, '-') << "\n";

    auto trunc = [](const std::string& s, int w) {
        return (int)s.size() <= w ? s : s.substr(0, w - 2) + "..";
    };

    for (std::size_t i = 0; i < solution.size(); i++) {
        const PrimerResult& r = solution[i];
        std::string pdr_str   = "[" + std::to_string(r.pdr_left)
                              + ", " + std::to_string(r.pdr_right) + "]";

        std::cout << std::right
                  << std::setw(W_IDX)  << i
                  << std::setw(W_POOL) << r.pool
                  << std::setw(W_SEG)  << r.segment_id
                  << std::setw(W_PDR)  << trunc(pdr_str, W_PDR)
                  << std::setw(W_SEQ)  << trunc(r.left.seq,  W_SEQ)
                  << std::setw(W_TM)   << std::fixed << std::setprecision(1) << r.left.tm
                  << std::setw(W_GC)   << std::fixed << std::setprecision(1) << r.left.gc
                  << std::setw(W_SEQ)  << trunc(r.right.seq, W_SEQ)
                  << std::setw(W_TM)   << std::fixed << std::setprecision(1) << r.right.tm
                  << std::setw(W_GC)   << std::fixed << std::setprecision(1) << r.right.gc
                  << std::setw(W_SIZE) << r.product_size
                  << "\n";
    }

    std::cout << std::string(TOTAL, '-') << "\n\n";
}

struct PoolDimerSummary {
    std::size_t amplicons = 0, affected_amplicons = 0, bad_pairs = 0;
    std::size_t affected_primers = 0;
    double severity = 0;
    double cutoff = 0;
};

static std::map<std::size_t, PoolDimerSummary> summarize_dimers(
    const std::vector<PrimerResult>& solution) {
    std::map<std::size_t, PoolDimerSummary> pools;
    for (const auto& pair : solution) {
        if (!pair.left.dimer_risk.cutoff || !pair.right.dimer_risk.cutoff) continue;
        auto& pool = pools[pair.pool];
        ++pool.amplicons;
        pool.affected_amplicons += !pair.left.dimer_risk.partners.empty() ||
                                  !pair.right.dimer_risk.partners.empty();
        pool.cutoff = *pair.left.dimer_risk.cutoff;
        for (const auto* primer : {&pair.left, &pair.right}) {
            pool.affected_primers += !primer->dimer_risk.partners.empty();
            pool.bad_pairs += primer->dimer_risk.partners.size();
            for (const auto& partner : primer->dimer_risk.partners)
                pool.severity += *primer->dimer_risk.cutoff - partner.dg;
        }
    }
    for (auto& entry : pools) {
        entry.second.bad_pairs /= 2;
        entry.second.severity /= 2; // Each pair appears at both primer endpoints.
    }
    return pools;
}

void annotate_solution_dimers(const std::vector<SegmentContext>& segments,
                              std::vector<PrimerResult>& solution,
                              double cutoff, double min_c, double max_c,
                              const std::string& output_prefix) {
    if (!std::isfinite(cutoff))
        throw std::invalid_argument("Dimer annotation cutoff must be finite");
    Thal::validate_dimer_range(min_c, max_c);
    std::vector<std::size_t> order(solution.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](std::size_t a, std::size_t b) {
        const auto& first = solution[a];
        const auto& second = solution[b];
        if (first.segment_id != second.segment_id)
            return first.segment_id < second.segment_id;
        return std::min(first.left.start, first.right.start) <
               std::min(second.left.start, second.right.start);
    });

    struct FinalPrimer {
        Oligo* oligo;
        std::string name;
        std::string sequence;
        std::size_t pool;
    };
    std::vector<FinalPrimer> primers(solution.size() * 2);
    std::map<std::size_t, std::size_t> segment_counts;
    for (const auto index : order) {
        auto& result = solution[index];
        const auto prefix = segments.at(result.segment_id).name + "_" +
                            std::to_string(++segment_counts[result.segment_id]);
        result.left.dimer_risk = {cutoff, {}};
        result.right.dimer_risk = {cutoff, {}};
        primers[2 * index] = {&result.left, prefix + "_LEFT_1", result.left.seq, result.pool};
        primers[2 * index + 1] = {&result.right, prefix + "_RIGHT_1", revcomp(result.right.seq), result.pool};
    }

    const auto path = output_prefix + ".bad-dimers.tsv";
    std::ofstream report(path);
    if (!report) throw std::runtime_error("Cannot write " + path);
    report << "pool\tprimer_a\tprimer_b\tdg_cal_per_mol\tcutoff_cal_per_mol"
           << "\ttemp_min_c\ttemp_max_c\n"
           << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (std::size_t i = 0; i < primers.size(); ++i) {
        for (std::size_t j = i + 1; j < primers.size(); ++j) {
            auto& a = primers[i];
            auto& b = primers[j];
            if (a.pool != b.pool) continue;
            const double dg = Thal::compute_dimer_worst_dg(
                a.sequence, b.sequence, min_c, max_c);
            if (!std::isfinite(dg))
                throw std::runtime_error("Non-finite selected primer dimer score");
            if (dg >= cutoff) continue;
            a.oligo->dimer_risk.partners.push_back({b.name, dg});
            b.oligo->dimer_risk.partners.push_back({a.name, dg});
            report << a.pool << '\t' << a.name << '\t' << b.name << '\t' << dg
                   << '\t' << cutoff << '\t' << min_c << '\t' << max_c << '\n';
        }
    }
    for (auto& primer : primers)
        std::sort(primer.oligo->dimer_risk.partners.begin(),
                  primer.oligo->dimer_risk.partners.end(),
                  [](const DimerPartner& a, const DimerPartner& b) { return a.dg < b.dg; });
    report.close();
    if (!report) throw std::runtime_error("Failed writing " + path);
    std::cout << "Selected below-cutoff dimer pairs written to: " << path << "\n";
    const auto summary_path = output_prefix + ".dimer-summary.tsv";
    std::ofstream summary(summary_path);
    if (!summary) throw std::runtime_error("Cannot write " + summary_path);
    summary << "pool\taffected_amplicons\tamplicons\taffected_amplicon_percent"
            << "\tseverity_cal_per_mol\tbelow_cutoff_pairs\tcutoff_cal_per_mol"
            << "\taffected_primers\tprimers\taffected_primer_percent\n"
            << std::setprecision(std::numeric_limits<double>::max_digits10);
    for (const auto& entry : summarize_dimers(solution)) {
        const auto& pool = entry.second;
        summary << entry.first << '\t' << pool.affected_amplicons << '\t'
                << pool.amplicons << '\t' << 100.0 * pool.affected_amplicons / pool.amplicons
                << '\t' << pool.severity << '\t' << pool.bad_pairs << '\t' << pool.cutoff
                << '\t' << pool.affected_primers << '\t' << 2 * pool.amplicons
                << '\t' << 100.0 * pool.affected_primers / (2 * pool.amplicons) << '\n';
    }
    summary.close();
    if (!summary) throw std::runtime_error("Failed writing " + summary_path);
    std::cout << "Dimer impact summary written to: " << summary_path << "\n";
}

void save_solution(std::vector<SegmentContext>& seg, std::vector<PrimerResult>& solution, const std::string& filename) {
    std::ofstream out(filename);
    if (!out) {
        std::cerr << "Error: could not open output file: " << filename << "\n";
        return;
    }

    // Sort by segment_id, then by the earliest start in the pair
    std::vector<PrimerResult> sorted = solution;
    std::sort(sorted.begin(), sorted.end(), [](const PrimerResult& a, const PrimerResult& b) {
        if (a.segment_id != b.segment_id)
            return a.segment_id < b.segment_id;
        std::size_t a_min = std::min(a.left.start, a.right.start);
        std::size_t b_min = std::min(b.left.start, b.right.start);
        return a_min < b_min;
    });

    // Flatten into individual primer entries with their strand
    struct PrimerEntry {
        std::size_t segment_id;
        std::size_t start;
        std::size_t end;
        std::size_t pool;
        std::string seq;
        char        strand;
        int         amp_idx;  // filled after per-seg counting
        PrimerDimerRisk dimer_risk;
    };

    std::vector<PrimerEntry> entries;
    std::unordered_map<std::size_t, int> seg_amp_count;

    for (const PrimerResult& r : sorted) {
        int& cnt = seg_amp_count[r.segment_id];
        cnt++;

        entries.push_back({ r.segment_id,
                            r.left.start, r.left.start + r.left.length,
                            r.pool, r.left.seq, '+', cnt, r.left.dimer_risk });

        entries.push_back({ r.segment_id,
                            r.right.start, r.right.start + r.right.length,
                            r.pool, r.right.seq, '-', cnt, r.right.dimer_risk });
    }

    // Sort all individual primers by segment then start position
    std::sort(entries.begin(), entries.end(), [](const PrimerEntry& a, const PrimerEntry& b) {
        if (a.segment_id != b.segment_id)
            return a.segment_id < b.segment_id;
        return a.start < b.start;
    });

    for (const PrimerEntry& e : entries) {
        std::string strand_name = (e.strand == '+') ? "LEFT" : "RIGHT";
        std::string name = seg[e.segment_id].name + "_"
                         + std::to_string(e.amp_idx)    + "_"
                         + strand_name                   + "_1";

        out << seg[e.segment_id].name << "\t"
            << e.start      << "\t"
            << e.end        << "\t"
            << name         << "\t"
            << e.pool       << "\t"
            << e.strand     << "\t"
            << e.seq        << "\t";
        if (!e.dimer_risk.cutoff) {
            out << "NOT_EVALUATED\t.\t.\t.\t.\n";
            continue;
        }
        const auto& partners = e.dimer_risk.partners;
        out << (partners.empty() ? "NONE_BELOW_CUTOFF" : "BELOW_CUTOFF") << '\t'
            << partners.size() << '\t'
            << std::setprecision(std::numeric_limits<double>::max_digits10);
        if (partners.empty()) out << '.';
        else out << partners.front().dg;
        out << '\t' << *e.dimer_risk.cutoff << '\t';
        if (partners.empty()) out << '.';
        for (std::size_t i = 0; i < partners.size(); ++i) {
            if (i) out << ',';
            out << partners[i].primer_name;
        }
        out << '\n';
    }

    std::cout << "Solution written to: " << filename << "\n";
}

static std::string xml_escape(const std::string& text) {
    std::string escaped;
    escaped.reserve(text.size());
    for (char c : text) {
        switch (c) {
            case '&':  escaped += "&amp;";  break;
            case '<':  escaped += "&lt;";   break;
            case '>':  escaped += "&gt;";   break;
            case '\"': escaped += "&quot;"; break;
            case '\'': escaped += "&apos;"; break;
            default:   escaped += c;        break;
        }
    }
    return escaped;
}

static std::string dimer_tooltip(const PrimerDimerRisk& risk) {
    if (!risk.cutoff) return "&#10;Dimer status: NOT_EVALUATED";
    std::ostringstream text;
    text << "&#10;Dimer status: "
         << (risk.partners.empty() ? "NONE_BELOW_CUTOFF" : "BELOW_CUTOFF")
         << "&#10;Dimer cutoff: " << std::fixed << std::setprecision(2)
         << *risk.cutoff << " cal/mol"
         << "&#10;Dimer partner count: " << risk.partners.size();
    if (!risk.partners.empty())
        text << "&#10;Dimer worst dG: " << risk.partners.front().dg << " cal/mol";
    const auto shown = std::min<std::size_t>(risk.partners.size(), 5);
    for (std::size_t i = 0; i < shown; ++i)
        text << "&#10;Dimer partner: " << xml_escape(risk.partners[i].primer_name)
             << " (" << risk.partners[i].dg << " cal/mol)";
    if (risk.partners.size() > shown)
        text << "&#10;Dimer partners remaining: " << risk.partners.size() - shown
             << " (see .bad-dimers.tsv)";
    return text.str();
}

void save_coverage_report(const std::vector<SegmentContext>& segments,
                          const std::vector<PrimerResult>& solution,
                          const std::string& output_prefix) {
    struct Amplicon {
        std::size_t index;
        std::size_t start;
        std::size_t end;
        std::size_t left_start;
        std::size_t left_end;
        std::size_t right_start;
        std::size_t right_end;
        std::size_t pool;
        Oligo left;
        Oligo right;
        bool valid;
        std::string problem;
    };

    const std::string tsv_name = output_prefix + ".coverage.tsv";
    const std::string svg_name = output_prefix + ".coverage.svg";
    const std::string variation_name = output_prefix + ".variation.tsv";
    std::ofstream tsv(tsv_name);
    std::ofstream svg(svg_name);
    std::ofstream variation_out(variation_name);
    if (!tsv || !svg || !variation_out) {
        std::cerr << "Error: could not create coverage report files\n";
        return;
    }

    constexpr double X0 = 170.0;
    constexpr double PIXELS_PER_BASE = 2.0;
    constexpr double RIGHT_MARGIN = 50.0;
    constexpr double ROW_HEIGHT = 105.0;
    constexpr double SEGMENT_HEADER_HEIGHT = 42.0;
    std::size_t maximum_template_length = 0;
    for (const auto& segment : segments) {
        maximum_template_length = std::max(
            maximum_template_length, segment.tmpl.size());
    }
    const double svg_width = X0 + RIGHT_MARGIN +
        std::max(1.0, maximum_template_length * PIXELS_PER_BASE);
    const std::size_t total_rows = segments.size();
    const auto pool_summaries = summarize_dimers(solution);
    const double summary_height = 16.0 * pool_summaries.size();
    const double svg_height = 61.0 + summary_height +
        total_rows * ROW_HEIGHT +
        segments.size() * (SEGMENT_HEADER_HEIGHT + 15.0);

    struct PrimerEndpoint {
        double x, y;
        std::size_t pool;
        PrimerDimerRisk risk;
    };
    std::map<std::string, PrimerEndpoint> primer_endpoints;

    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\""
        << svg_width << "\" height=\"" << svg_height
        << "\" viewBox=\"0 0 " << svg_width << " " << svg_height << "\">\n"
        << "<style>text{font-family:Arial,sans-serif;font-size:12px}"
        << ".title{font-size:18px;font-weight:bold}.label{font-weight:bold}"
        << ".axis{stroke:#334155;stroke-width:1}.gap{fill:#ef4444;opacity:.75}"
        << ".p1{fill:#2563eb;stroke:#1d4ed8}.p2{fill:#16a34a;stroke:#15803d}"
        << ".primer{cursor:pointer}.primer-label{font-size:10px;pointer-events:none}"
        << ".bad-dimer{stroke:#f97316;stroke-width:2.5}"
        << ".bad-dimer-label{fill:#c2410c;font-weight:bold}"
        << ".dimer-link{opacity:.4}.dimer-connection{fill:none;stroke:#ea580c;stroke-width:1.8;pointer-events:stroke}"
        << ".dimer-endpoint{fill:#ea580c;pointer-events:none}"
        << ".dimer-link:hover{opacity:1}.dimer-link:hover .dimer-connection{stroke-width:3.5}"
        << ".dimer-focus .dimer-link{opacity:.04}.dimer-focus .dimer-link.active{opacity:1}"
        << ".dimer-focus .dimer-link.active .dimer-connection{stroke-width:3}"
        << ".primer:hover{opacity:.78}.primer-tooltip text{font-size:12px;fill:#f8fafc}"
        << ".primer-tooltip .heading{font-weight:bold;font-size:13px}"
        << ".primer-tooltip .sequence{font-family:Consolas,monospace;white-space:pre}"
        << ".invalid{fill:#9ca3af;stroke:#dc2626;stroke-width:2}</style>\n"
        << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n"
        << "<text x=\"20\" y=\"27\" class=\"title\">Multiplex amplicon coverage</text>\n"
        << "<text x=\"430\" y=\"27\">Variation: pale = conserved, red = variable, gray = no valid calls</text>\n"
        << "<text x=\"20\" y=\"45\" fill=\"#c2410c\">Orange curves connect below-cutoff dimer pairs; * marks affected primers (hover to highlight)</text>\n";
    double summary_y = 63.0;
    for (const auto& entry : pool_summaries) {
        const auto& pool = entry.second;
        svg << "<text x=\"20\" y=\"" << summary_y << "\">Pool "
            << entry.first << " — affected primers " << pool.affected_primers << "/"
            << 2 * pool.amplicons << " (" << std::fixed << std::setprecision(1)
            << 100.0 * pool.affected_primers / (2 * pool.amplicons)
            << "%); amplicons "
            << pool.affected_amplicons << "/" << pool.amplicons << " ("
            << std::fixed << std::setprecision(1)
            << 100.0 * pool.affected_amplicons / pool.amplicons
            << "%); severity " << std::setprecision(2) << pool.severity / 1000.0
            << " kcal/mol</text>\n";
        summary_y += 16.0;
    }

    tsv << "record\tsegment\tindex\tstart\tend\tlength\tpool\tvalid"
        << "\tcovered_bases\tcoverage_percent\tgap_count\tproblem\n";
    variation_out
        << "segment\tposition\ttemplate_base\tvalid_sequences\ttotal_sequences"
        << "\tvariation_rate\tmissing_rate\tA\tC\tG\tT\n";

    double segment_top = 64.0 + summary_height;
    for (std::size_t segment_id = 0; segment_id < segments.size(); ++segment_id) {
        const auto& segment = segments[segment_id];
        const std::size_t template_length = segment.tmpl.size();
        std::vector<Amplicon> amplicons;

        struct SiteVariation {
            std::size_t valid = 0;
            std::size_t counts[4] = {0, 0, 0, 0};
            double variation = 0.0;
            double missing = 1.0;
        };
        std::vector<SiteVariation> site_variation(template_length);

        auto base_index = [](char base) -> int {
            switch (static_cast<char>(std::toupper(
                        static_cast<unsigned char>(base)))) {
                case 'A': return 0;
                case 'C': return 1;
                case 'G': return 2;
                case 'T': return 3;
                default:  return -1;
            }
        };

        for (std::size_t pos = 0; pos < template_length; ++pos) {
            SiteVariation& site = site_variation[pos];
            for (const std::string& sequence : segment.sequences) {
                if (pos >= sequence.size()) continue;
                const int index = base_index(sequence[pos]);
                if (index < 0) continue;
                ++site.counts[index];
                ++site.valid;
            }

            if (site.valid > 0) {
                const std::size_t major = *std::max_element(
                    std::begin(site.counts), std::end(site.counts));
                site.variation = 1.0 -
                    static_cast<double>(major) /
                    static_cast<double>(site.valid);
            }
            if (!segment.sequences.empty()) {
                site.missing = 1.0 -
                    static_cast<double>(site.valid) /
                    static_cast<double>(segment.sequences.size());
            }

            variation_out << segment.name << "\t" << pos << "\t"
                << (pos < segment.tmpl.size() ? segment.tmpl[pos] : 'N') << "\t"
                << site.valid << "\t" << segment.sequences.size() << "\t"
                << std::fixed << std::setprecision(6) << site.variation << "\t"
                << site.missing << "\t" << site.counts[0] << "\t"
                << site.counts[1] << "\t" << site.counts[2] << "\t"
                << site.counts[3] << "\n";
        }

        for (const PrimerResult& result : solution) {
            if (result.segment_id != segment_id) continue;

            Amplicon amp{};
            amp.left_start  = result.left.start;
            amp.left_end    = result.left.start + result.left.length;
            amp.right_start = result.right.start;
            amp.right_end   = result.right.start + result.right.length;
            amp.start       = amp.left_start;
            amp.end         = amp.right_end;
            amp.pool        = result.pool;
            amp.left        = result.left;
            amp.right       = result.right;
            amp.valid       = true;

            if (result.left.length == 0 || result.right.length == 0) {
                amp.valid = false;
                amp.problem = "zero-length primer";
            } else if (amp.left_start >= amp.right_start) {
                amp.valid = false;
                amp.problem = "invalid primer orientation/order";
            } else if (amp.left_end > template_length ||
                       amp.right_end > template_length) {
                amp.valid = false;
                amp.problem = "primer outside template bounds";
            }
            amplicons.push_back(std::move(amp));
        }

        std::sort(amplicons.begin(), amplicons.end(),
                  [](const Amplicon& a, const Amplicon& b) {
                      if (a.start != b.start) return a.start < b.start;
                      return a.end < b.end;
                  });

        // Match save_solution() BED numbering: genomic order within segment.
        for (std::size_t i = 0; i < amplicons.size(); ++i) {
            amplicons[i].index = i + 1;
        }

        std::vector<std::pair<std::size_t, std::size_t>> merged;
        for (const Amplicon& amp : amplicons) {
            if (!amp.valid) continue;
            if (merged.empty() || amp.start > merged.back().second) {
                merged.push_back({amp.start, amp.end});
            } else {
                merged.back().second = std::max(merged.back().second, amp.end);
            }
        }

        std::vector<std::pair<std::size_t, std::size_t>> gaps;
        std::size_t cursor = 0;
        std::size_t covered = 0;
        for (const auto& interval : merged) {
            if (cursor < interval.first) gaps.push_back({cursor, interval.first});
            covered += interval.second - interval.first;
            cursor = interval.second;
        }
        if (cursor < template_length) gaps.push_back({cursor, template_length});

        const double coverage_percent = template_length == 0
            ? 0.0
            : 100.0 * static_cast<double>(covered) /
              static_cast<double>(template_length);

        tsv << "SUMMARY\t" << segment.name << "\t.\t0\t"
            << template_length << "\t" << template_length << "\t.\t.\t"
            << covered << "\t" << std::fixed << std::setprecision(3)
            << coverage_percent << "\t" << gaps.size() << "\t.\n";

        for (const Amplicon& amp : amplicons) {
            tsv << "AMPLICON\t" << segment.name << "\t" << amp.index << "\t"
                << amp.start << "\t" << amp.end << "\t"
                << (amp.end >= amp.start ? amp.end - amp.start : 0) << "\t"
                << amp.pool << "\t" << (amp.valid ? 1 : 0)
                << "\t.\t.\t.\t" << (amp.problem.empty() ? "." : amp.problem)
                << "\n";
        }
        for (const auto& gap : gaps) {
            tsv << "GAP\t" << segment.name << "\t.\t" << gap.first << "\t"
                << gap.second << "\t" << gap.second - gap.first
                << "\t.\t.\t.\t.\t.\tuncovered\n";
        }

        svg << "<text x=\"20\" y=\"" << segment_top + 14
            << "\" class=\"label\">" << xml_escape(segment.name) << "</text>\n"
            << "<text x=\"20\" y=\"" << segment_top + 31
            << "\">" << covered << "/" << template_length << " bp ("
            << std::fixed << std::setprecision(2) << coverage_percent
            << "%), gaps: " << gaps.size() << "</text>\n";

        // Draw the entire template on one continuous axis. The SVG canvas is
        // intentionally wider than the viewport so it can be inspected with
        // horizontal scrolling at per-primer resolution.
        constexpr std::size_t row_count = 1;
        const double plot_width =
            std::max(1.0, template_length * PIXELS_PER_BASE);

        for (std::size_t row = 0; row < row_count; ++row) {
            const std::size_t row_start = 0;
            const std::size_t row_end = template_length;
            const std::size_t row_length = std::max<std::size_t>(1, row_end - row_start);
            const double base_y = segment_top + SEGMENT_HEADER_HEIGHT + row * ROW_HEIGHT;

            auto xpos = [&](std::size_t position) {
                const std::size_t clipped = std::min(
                    row_end, std::max(row_start, position));
                return X0 + plot_width *
                    static_cast<double>(clipped - row_start) /
                    static_cast<double>(row_length);
            };

            svg << "<line x1=\"" << X0 << "\" y1=\"" << base_y
                << "\" x2=\"" << X0 + plot_width << "\" y2=\""
                << base_y << "\" class=\"axis\"/>\n";

            // Per-site variation heatmap. The theoretical maximum for four
            // equally frequent bases is 0.75, which maps to saturated red.
            for (std::size_t pos = row_start; pos < row_end; ++pos) {
                const SiteVariation& site = site_variation[pos];
                const double normalized = std::min(1.0, site.variation / 0.75);
                const int red   = 255;
                const int green = static_cast<int>(255.0 * (1.0 - normalized));
                const int blue  = static_cast<int>(220.0 * (1.0 - normalized));
                const double x1 = xpos(pos);
                const double x2 = xpos(pos + 1);

                svg << "<rect x=\"" << x1 << "\" y=\"" << base_y + 10
                    << "\" width=\"" << std::max(0.05, x2 - x1)
                    << "\" height=\"10\" fill=\"";
                if (site.valid == 0) svg << "#94a3b8";
                else svg << "rgb(" << red << "," << green << "," << blue << ")";
                svg << "\"><title>Position " << pos
                    << "; template " << xml_escape(std::string(1, segment.tmpl[pos]))
                    << "; variation " << std::fixed << std::setprecision(4)
                    << site.variation << "; missing " << site.missing
                    << "; A/C/G/T " << site.counts[0] << "/" << site.counts[1]
                    << "/" << site.counts[2] << "/" << site.counts[3]
                    << "</title></rect>\n";
            }

            for (const auto& gap : gaps) {
                const std::size_t start = std::max(gap.first, row_start);
                const std::size_t end = std::min(gap.second, row_end);
                if (start >= end) continue;
                svg << "<rect x=\"" << xpos(start) << "\" y=\""
                    << base_y - 5 << "\" width=\""
                    << std::max(1.0, xpos(end) - xpos(start))
                    << "\" height=\"10\" class=\"gap\"><title>Gap "
                    << gap.first << "-" << gap.second << " ("
                    << gap.second - gap.first << " bp)</title></rect>\n";
            }

            for (const Amplicon& amp : amplicons) {
                const std::size_t visible_start = std::max(amp.start, row_start);
                const std::size_t visible_end = std::min(amp.end, row_end);
                if (visible_start >= visible_end) continue;

                const double y = base_y + (amp.pool == 2 ? 68.0 : 41.0);
                const std::string cls = amp.valid
                    ? (amp.pool == 2 ? "p2" : "p1") : "invalid";
                svg << "<rect x=\"" << xpos(visible_start) << "\" y=\"" << y
                    << "\" width=\""
                    << std::max(1.0, xpos(visible_end) - xpos(visible_start))
                    << "\" height=\"12\" rx=\"2\" class=\"" << cls
                    << "\"><title>Amplicon " << amp.index << ", pool " << amp.pool
                    << ": " << amp.start << "-" << amp.end;
                if (!amp.problem.empty()) svg << "; " << xml_escape(amp.problem);
                svg << "</title></rect>\n";

                const auto name_prefix = segment.name + "_" + std::to_string(amp.index);
                auto primer_detail_attributes = [&](const Oligo& oligo,
                                                     std::size_t start,
                                                     std::size_t end,
                                                     bool is_right) {
                    const std::string primer_sequence = is_right
                        ? revcomp(oligo.seq) : oligo.seq;
                    std::string target_sequence;
                    target_sequence.reserve(primer_sequence.size());
                    for (char base : primer_sequence) {
                        target_sequence.push_back(comp(base));
                    }

                    std::ostringstream variation;
                    std::ostringstream positions;
                    variation << std::fixed << std::setprecision(1);
                    for (std::size_t j = 0; j < primer_sequence.size(); ++j) {
                        if (j != 0) {
                            variation << ',';
                            positions << ',';
                        }
                        const std::size_t position = is_right
                            ? end - 1 - j : start + j;
                        positions << position;
                        variation << (position < site_variation.size()
                            ? 100.0 * site_variation[position].variation : 0.0);
                    }

                    std::ostringstream attributes;
                    attributes << " data-primer-id=\""
                               << xml_escape(name_prefix + (is_right ? "_RIGHT_1" : "_LEFT_1"))
                               << "\" data-primer-seq=\""
                               << xml_escape(primer_sequence)
                               << "\" data-target-seq=\""
                               << xml_escape(target_sequence)
                               << "\" data-site-variation=\""
                               << variation.str()
                               << "\" data-site-positions=\""
                               << positions.str() << "\"";
                    return attributes.str();
                };

                if (amp.valid && amp.left_end > row_start && amp.left_start < row_end) {
                    const double lx1 = xpos(amp.left_start);
                    const double lx2 = xpos(amp.left_end);
                    primer_endpoints[name_prefix + "_LEFT_1"] = {
                        (lx1 + lx2) / 2.0, y + 6, amp.pool, amp.left.dimer_risk};
                    svg << "<polygon points=\"" << lx1 << "," << y - 3 << " "
                        << lx1 << "," << y + 15 << " " << lx2 << "," << y + 6
                        << "\" fill=\"#0f172a\" class=\"primer"
                        << (amp.left.dimer_risk.partners.empty() ? "" : " bad-dimer") << "\""
                        << primer_detail_attributes(amp.left, amp.left_start,
                                                    amp.left_end, false)
                        << "><title>LEFT primer "
                        << amp.index << "&#10;Pool: " << amp.pool
                        << "&#10;Sequence: " << xml_escape(amp.left.seq)
                        << "&#10;Coordinates: " << amp.left_start << "-"
                        << amp.left_end << "&#10;Length: " << amp.left.length
                        << "&#10;Tm: " << std::fixed << std::setprecision(2)
                        << amp.left.tm << " C&#10;GC: " << amp.left.gc
                        << "%"
                        << dimer_tooltip(amp.left.dimer_risk)
                        << "</title></polygon>\n";
                    svg << "<text x=\"" << (lx1 + lx2) / 2.0
                        << "\" y=\"" << y - 5
                        << "\" text-anchor=\"middle\" class=\"primer-label"
                        << (amp.left.dimer_risk.partners.empty() ? "" : " bad-dimer-label") << "\">L"
                        << amp.index << (amp.left.dimer_risk.partners.empty() ? "" : "*") << "</text>\n";
                }
                if (amp.valid && amp.right_end > row_start && amp.right_start < row_end) {
                    const double rx1 = xpos(amp.right_start);
                    const double rx2 = xpos(amp.right_end);
                    primer_endpoints[name_prefix + "_RIGHT_1"] = {
                        (rx1 + rx2) / 2.0, y + 6, amp.pool, amp.right.dimer_risk};
                    svg << "<polygon points=\"" << rx2 << "," << y - 3 << " "
                        << rx2 << "," << y + 15 << " " << rx1 << "," << y + 6
                        << "\" fill=\"#0f172a\" class=\"primer"
                        << (amp.right.dimer_risk.partners.empty() ? "" : " bad-dimer") << "\""
                        << primer_detail_attributes(amp.right, amp.right_start,
                                                    amp.right_end, true)
                        << "><title>RIGHT primer "
                        << amp.index << "&#10;Pool: " << amp.pool
                        << "&#10;Sequence: " << xml_escape(amp.right.seq)
                        << "&#10;Coordinates: " << amp.right_start << "-"
                        << amp.right_end << "&#10;Length: " << amp.right.length
                        << "&#10;Tm: " << std::fixed << std::setprecision(2)
                        << amp.right.tm << " C&#10;GC: " << amp.right.gc
                        << "%"
                        << dimer_tooltip(amp.right.dimer_risk)
                        << "</title></polygon>\n";
                    svg << "<text x=\"" << (rx1 + rx2) / 2.0
                        << "\" y=\"" << y - 5
                        << "\" text-anchor=\"middle\" class=\"primer-label"
                        << (amp.right.dimer_risk.partners.empty() ? "" : " bad-dimer-label") << "\">R"
                        << amp.index << (amp.right.dimer_risk.partners.empty() ? "" : "*") << "</text>\n";
                }
            }

            svg << "<text x=\"105\" y=\"" << base_y + 19
                << "\">Variation</text>\n"
                << "<text x=\"120\" y=\"" << base_y + 52
                << "\" fill=\"#2563eb\">Pool 1</text>\n"
                << "<text x=\"120\" y=\"" << base_y + 79
                << "\" fill=\"#16a34a\">Pool 2</text>\n"
                << "<text x=\"" << X0 << "\" y=\"" << base_y + 100
                << "\">" << row_start << "</text><text x=\""
                << X0 + plot_width << "\" y=\"" << base_y + 100
                << "\" text-anchor=\"end\">" << row_end << " bp</text>\n";
        }

        segment_top += SEGMENT_HEADER_HEIGHT + row_count * ROW_HEIGHT + 15.0;
    }

    svg << "<g id=\"dimer-connections\">\n";
    for (const auto& entry : primer_endpoints) {
        const auto& name = entry.first;
        const auto& first = entry.second;
        for (const auto& partner : first.risk.partners) {
            if (name >= partner.primer_name) continue;
            const auto target = primer_endpoints.find(partner.primer_name);
            if (target == primer_endpoints.end() || first.pool != target->second.pool) continue;
            const auto& second = target->second;
            const double control_y = std::min(first.y, second.y) - 32.0;
            svg << "<g class=\"dimer-link\" data-primer-a=\"" << xml_escape(name)
                << "\" data-primer-b=\"" << xml_escape(partner.primer_name)
                << "\" data-dg=\"" << std::setprecision(8) << partner.dg << "\">"
                << "<title>" << xml_escape(name) << " &lt;-&gt; " << xml_escape(partner.primer_name)
                << "&#10;Pool: " << first.pool << "&#10;Dimer dG: "
                << std::fixed << std::setprecision(2) << partner.dg << " cal/mol";
            if (first.risk.cutoff)
                svg << "&#10;Cutoff: " << *first.risk.cutoff << " cal/mol";
            svg << "</title><path class=\"dimer-connection\" d=\"M "
                << first.x << " " << first.y << " C " << first.x << " " << control_y
                << ", " << second.x << " " << control_y << ", "
                << second.x << " " << second.y << "\"/>"
                << "<circle class=\"dimer-endpoint\" cx=\"" << first.x
                << "\" cy=\"" << first.y << "\" r=\"2.5\"/>"
                << "<circle class=\"dimer-endpoint\" cx=\"" << second.x
                << "\" cy=\"" << second.y << "\" r=\"2.5\"/></g>\n";
        }
    }
    svg << "</g>\n";

    // A native SVG <title> is useful as a no-script fallback, but its tooltip
    // is delayed and varies by viewer.  Upgrade primer titles to a consistent,
    // pointer-following tooltip when the SVG is opened in a browser.
    svg << R"SVG(<g id="primer-tooltip" class="primer-tooltip" visibility="hidden" pointer-events="none">
<rect x="0" y="0" width="10" height="10" rx="4" fill="#0f172a" opacity="0.96"/>
<text x="10" y="18"></text>
</g>
<script><![CDATA[
(function () {
  const root = document.documentElement;
  const tooltip = document.getElementById('primer-tooltip');
  const box = tooltip.querySelector('rect');
  const text = tooltip.querySelector('text');
  const padding = 10;
  const lineHeight = 17;
  const dimerLinks = Array.from(root.querySelectorAll('.dimer-link'));

  function point(event) {
    const p = root.createSVGPoint();
    p.x = event.clientX;
    p.y = event.clientY;
    return p.matrixTransform(root.getScreenCTM().inverse());
  }

  function visibleBounds() {
    const rect = root.getBoundingClientRect();
    const matrix = root.getScreenCTM().inverse();
    const corner = function (clientX, clientY) {
      const p = root.createSVGPoint();
      p.x = clientX;
      p.y = clientY;
      return p.matrixTransform(matrix);
    };
    const topLeft = corner(Math.max(0, rect.left) + 4,
                           Math.max(0, rect.top) + 4);
    const bottomRight = corner(Math.min(window.innerWidth, rect.right) - 4,
                               Math.min(window.innerHeight, rect.bottom) - 4);
    return {
      left: topLeft.x,
      top: topLeft.y,
      right: bottomRight.x,
      bottom: bottomRight.y
    };
  }

  function move(event) {
    const p = point(event);
    const bounds = text.getBBox();
    const visible = visibleBounds();
    const width = bounds.width + padding * 2;
    const height = bounds.height + padding * 2;
    let x = p.x + 14;
    let y = p.y + 14;
    if (x + width > visible.right) x = p.x - width - 14;
    if (y + height > visible.bottom) y = p.y - height - 14;
    x = Math.max(visible.left, Math.min(x, visible.right - width));
    y = Math.max(visible.top, Math.min(y, visible.bottom - height));
    tooltip.setAttribute('transform', `translate(${x},${y})`);
  }

  document.querySelectorAll('.primer').forEach(function (primer) {
    const title = primer.querySelector('title');
    if (!title) return;
    const detail = title.textContent;
    primer.setAttribute('aria-label', detail.replace(/\n/g, '; '));
    title.remove();

    primer.addEventListener('pointerenter', function (event) {
      root.classList.add('dimer-focus');
      dimerLinks.forEach(function (link) {
        link.classList.toggle('active',
          link.dataset.primerA === primer.dataset.primerId ||
          link.dataset.primerB === primer.dataset.primerId);
      });
      text.replaceChildren();
      const metadata = Object.fromEntries(detail.split('\n').slice(1).map(function (line) {
        const separator = line.indexOf(':');
        return separator < 0 ? [line, ''] : [line.slice(0, separator), line.slice(separator + 1).trim()];
      }));
      const primerBases = (primer.dataset.primerSeq || '').split('');
      const targetBases = (primer.dataset.targetSeq || '').split('');
      const positions = (primer.dataset.sitePositions || '').split(',');
      const variation = (primer.dataset.siteVariation || '').split(',');
      const cellWidth = 1 + Math.max(3, ...positions.concat(variation).map(function (value) {
        return String(value).length;
      }));
      const availableWidth = Math.max(240, window.innerWidth - 32);
      const columnsPerBlock = Math.max(4, Math.floor(
        (availableWidth - 100) / (cellWidth * 7.25)));
      const tokenRow = function (label, values) {
        return label.padEnd(12) + values.map(function (value) {
          return String(value).padStart(cellWidth);
        }).join('');
      };
      const lines = [
        detail.split('\n')[0],
        'Pool ' + metadata.Pool + '  |  Coordinates ' + metadata.Coordinates,
        'Tm ' + metadata.Tm + '  |  GC ' + metadata.GC
      ];
      lines.push(...detail.split('\n').filter(function (line) {
        return line.startsWith('Dimer ');
      }));
      for (let start = 0; start < primerBases.length; start += columnsPerBlock) {
        const end = start + columnsPerBlock;
        lines.push(tokenRow('Position', positions.slice(start, end)));
        lines.push(tokenRow("Primer 5'", primerBases.slice(start, end)));
        lines.push(tokenRow("Target 3'", targetBases.slice(start, end)));
        lines.push(tokenRow('Variation %', variation.slice(start, end)));
      }
      lines.forEach(function (line, index) {
        const tspan = document.createElementNS('http://www.w3.org/2000/svg', 'tspan');
        tspan.setAttribute('x', padding);
        tspan.setAttribute('y', padding + 13 + index * lineHeight);
        if (index === 0) tspan.setAttribute('class', 'heading');
        if (index >= 2) {
          tspan.setAttribute('class', 'sequence');
          tspan.setAttributeNS('http://www.w3.org/XML/1998/namespace',
                               'xml:space', 'preserve');
        }
        tspan.textContent = line;
        text.appendChild(tspan);
      });
      const bounds = text.getBBox();
      box.setAttribute('width', bounds.width + padding * 2);
      box.setAttribute('height', bounds.height + padding * 2);
      tooltip.setAttribute('visibility', 'visible');
      move(event);
    });
    primer.addEventListener('pointermove', move);
    primer.addEventListener('pointerleave', function () {
      root.classList.remove('dimer-focus');
      tooltip.setAttribute('visibility', 'hidden');
    });
  });
}());
]]></script>
</svg>
)SVG";
    std::cout << "Coverage plot written to: " << svg_name << "\n"
              << "Coverage audit written to: " << tsv_name << "\n"
              << "Site variation written to: " << variation_name << "\n";
}
