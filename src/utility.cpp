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

void test_primer3() {
    // --- 1. Global settings (primer constraints) ---
    p3_global_settings *pa = p3_create_global_settings();
    p3_set_gs_primer_opt_size(pa, 20);
    p3_set_gs_primer_min_size(pa, 18);
    p3_set_gs_primer_max_size(pa, 25);
    p3_set_gs_primer_opt_tm(pa, 60.0);
    p3_set_gs_primer_min_tm(pa, 57.0);
    p3_set_gs_primer_max_tm(pa, 63.0);
    p3_set_gs_primer_min_gc(pa, 40.0);
    p3_set_gs_primer_max_gc(pa, 60.0);

    pa->p_args.salt_conc     = 50.0;
    pa->p_args.divalent_conc = 1.5;
    pa->p_args.dntp_conc     = 0.6;
    pa->p_args.dna_conc      = 50.0;

    // Product size range: 100–300 bp
    pa->pr_min[0] = 100;
    pa->pr_max[0] = 300;
    pa->num_intervals = 1;

    // Pick left + right primers (standard PCR pair)
    pa->pick_left_primer  = 1;
    pa->pick_right_primer = 1;
    pa->pick_internal_oligo = 0;

    // --- 2. Sequence arguments ---
    seq_args *sa = create_seq_arg();

    const char *tmpl = "GCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCC"
                       "TACATTTTAGCATCAGTGAGTACAGCATGCTTACTGG"
                       "AAGAGAGGGTCATGCAACAGATTAGGAGGTAAGTTTG"
                       "CAAAGGCAGGCTAAGGAGGAGACGCACTGAATGCCAT"
                       "GGTAAGAACTCTGGACATAAAAATATTGGAAGTTGTTG";

    sa->sequence = strdup(tmpl);
    sa->sequence_name = strdup("my_template");

    // Included region: start=0, length=full sequence
    sa->incl_s = 0;
    sa->incl_l = strlen(tmpl);

    // Number of primer pairs to return
    pa->num_return = 5;

    // --- 3. Run primer design ---
    p3retval *retval = choose_primers(pa, sa);

    if (retval == nullptr) {
        std::cerr << "choose_primers() returned NULL\n";
        return ;
    }

    // Check for errors
    if (retval->glob_err.data) {
        std::cerr << "Global error: " << retval->glob_err.data << "\n";
    }
    if (retval->per_sequence_err.data) {
        std::cerr << "Sequence error: " << retval->per_sequence_err.data << "\n";
    }

    // --- 4. Extract results ---
    int n_pairs = retval->best_pairs.num_pairs;
    std::cout << "Found " << n_pairs << " primer pair(s)\n\n";

    for (int i = 0; i < n_pairs; i++) {
        primer_pair *pp = &retval->best_pairs.pairs[i];
        primer_rec  *left  = pp->left;
        primer_rec  *right = pp->right;

        // Extract left primer sequence from template
        char left_seq[64] = {0};
        strncpy(left_seq, sa->sequence + left->start, left->length);

        // Right primer is on reverse complement — extract raw coords
        char right_seq[64] = {0};
        // right->start is the 3' end position; length goes leftward
        int right_start = right->start - right->length + 1;
        strncpy(right_seq, sa->sequence + right_start, right->length);

        std::cout << "Pair " << i + 1 << ":\n"
                  << "  Left:  " << left_seq
                  << "  (pos=" << left->start
                  << " len=" << left->length
                  << " Tm=" << left->temp << ")\n"
                  << "  Right: " << right_seq
                  << "  (pos=" << right->start
                  << " len=" << right->length
                  << " Tm=" << right->temp << ")\n"
                  << "  Product size: " << pp->product_size << "\n\n";
    }

    // --- 5. Cleanup ---
    destroy_p3retval(retval);
    destroy_seq_args(sa);
    p3_destroy_global_settings(pa);
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
