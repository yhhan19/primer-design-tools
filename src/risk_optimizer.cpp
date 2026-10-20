#include "risk_optimizer.hpp"

RiskOptimizer::RiskOptimizer(std::vector<std::string> &labels, 
                             std::vector<std::string> &data,
                             index_t len,
                             index_t min,
                             index_t max,
                             index_t fumax,
                             index_t rumax,
                             const Args &args) {
    this->args = args;
    std::vector<risk_t> input = compute_risk(data);
    this->tmpl = msa_consensus(data);
    size = input.size();
    risk = new risk_t[size];
    prefix_sum = new risk_t[size + 1];
    prefix_sum[0] = 0;
    for (index_t i = 0; i < input.size(); i ++) {
        risk[i] = input[i];
        prefix_sum[i + 1] = prefix_sum[i] + risk[i];
    }
    this->len = len;
    this->min = min;
    this->max = max;
    this->fumax = fumax;
    this->rumax = rumax;
    memo = new risk_t[KEY_LIMIT];
    prev = new key_t[KEY_LIMIT];
    solution_vector = new std::vector<std::pair<index_t, risk_t>>[size];
    solutions = new BST*[size];
    build_primer3_candidate_counts();
}

RiskOptimizer::RiskOptimizer(std::vector<risk_t> input, 
                             index_t len, 
                             index_t min, 
                             index_t max) {
    size = input.size();
    risk = new risk_t[size];
    prefix_sum = new risk_t[size + 1];
    prefix_sum[0] = 0;
    for (index_t i = 0; i < input.size(); i ++) {
        risk[i] = input[i];
        prefix_sum[i + 1] = prefix_sum[i] + risk[i];
    }
    this->len = len;
    this->min = min;
    this->max = max;
    memo = new risk_t[KEY_LIMIT];
    prev = new key_t[KEY_LIMIT];
    solution_vector = new std::vector<std::pair<index_t, risk_t>>[size];
    solutions = new BST*[size];
}

risk_t RiskOptimizer::cost(index_t p, risk_t u, risk_t alpha) {
    risk_t sum = prefix_sum[p + len] - prefix_sum[p];
    risk_t sum2 = sum * sum - u;
    if (sum2 < 0) return u * alpha;
    return sum2 + u * alpha;
}

void RiskOptimizer::build_primer3_candidate_counts() {
    std::cout << "Enumerating exact Primer3 candidates for all PDRs...\n";
    if (tmpl.size() < len) {
        throw std::runtime_error("template shorter than PDR length");
    }

    p3_global_settings* pa = p3_create_global_settings();
    if (pa == nullptr) {
        throw std::runtime_error("could not create Primer3 global settings");
    }

    p3_set_gs_primer_opt_size(pa, args.p3_opt_size);
    p3_set_gs_primer_min_size(pa, args.p3_min_size);
    p3_set_gs_primer_max_size(pa, args.p3_max_size);
    p3_set_gs_primer_opt_tm(pa, args.p3_opt_tm);
    p3_set_gs_primer_min_tm(pa, args.p3_min_tm);
    p3_set_gs_primer_max_tm(pa, args.p3_max_tm);
    p3_set_gs_primer_min_gc(pa, args.p3_min_gc);
    p3_set_gs_primer_max_gc(pa, args.p3_max_gc);
    pa->p_args.salt_conc     = args.mv;
    pa->p_args.divalent_conc = args.dv;
    pa->p_args.dntp_conc     = args.dntp;
    pa->p_args.dna_conc      = args.dna_conc;
    pa->primer_task          = pick_primer_list;
    pa->pick_left_primer     = 1;
    pa->pick_right_primer    = 1;
    pa->pick_internal_oligo  = 0;

    seq_args* sa = create_seq_arg();
    sa->sequence = strdup(tmpl.c_str());
    sa->sequence_name = strdup("pdr_candidate_preprocessing");

    p3retval* retval = choose_primers(pa, sa);
    if (retval == nullptr) {
        destroy_seq_args(sa);
        p3_destroy_global_settings(pa);
        throw std::runtime_error("Primer3 returned null during PDR preprocessing");
    }

    std::string error;
    if (retval->glob_err.data) error += retval->glob_err.data;
    if (retval->per_sequence_err.data) {
        if (!error.empty()) error += "; ";
        error += retval->per_sequence_err.data;
    }
    if (!error.empty()) {
        destroy_p3retval(retval);
        destroy_seq_args(sa);
        p3_destroy_global_settings(pa);
        throw std::runtime_error("Primer3 PDR preprocessing error: " + error);
    }
    if (retval->warnings.data) {
        std::cerr << "Primer3 PDR preprocessing warning: "
                  << retval->warnings.data << "\n";
    }

    const std::size_t num_pdrs = tmpl.size() - len + 1;
    std::vector<long long> left_diff(num_pdrs + 1, 0);
    std::vector<long long> right_diff(num_pdrs + 1, 0);

    auto add_to_containing_pdrs = [&](std::vector<long long>& diff,
                                      std::size_t primer_start,
                                      std::size_t primer_end) {
        if (primer_end <= primer_start || primer_end - primer_start > len) return;
        const std::size_t first = primer_end > len ? primer_end - len : 0;
        const std::size_t last = std::min(primer_start, num_pdrs - 1);
        if (first > last) return;
        ++diff[first];
        --diff[last + 1];
    };

    for (int i = 0; i < retval->fwd.num_elem; ++i) {
        const primer_rec& primer = retval->fwd.oligo[i];
        if (primer.start < 0 || primer.length <= 0) continue;
        const std::size_t start = static_cast<std::size_t>(primer.start);
        add_to_containing_pdrs(left_diff, start, start + primer.length);
    }
    for (int i = 0; i < retval->rev.num_elem; ++i) {
        const primer_rec& primer = retval->rev.oligo[i];
        if (primer.start < 0 || primer.length <= 0 ||
            primer.start - primer.length + 1 < 0) continue;
        const std::size_t start = static_cast<std::size_t>(
            primer.start - primer.length + 1);
        add_to_containing_pdrs(right_diff, start, start + primer.length);
    }

    left_valid_results.resize(num_pdrs);
    right_valid_results.resize(num_pdrs);
    long long left_running = 0;
    long long right_running = 0;
    for (std::size_t p = 0; p < num_pdrs; ++p) {
        left_running += left_diff[p];
        right_running += right_diff[p];
        left_valid_results[p] = static_cast<index_t>(left_running);
        right_valid_results[p] = static_cast<index_t>(right_running);
    }

    auto adaptive_minimum = [](const std::vector<index_t>& counts) {
        std::vector<index_t> nonzero;
        nonzero.reserve(counts.size());
        std::copy_if(counts.begin(), counts.end(),
                     std::back_inserter(nonzero),
                     [](index_t count) { return count > 0; });
        if (nonzero.empty()) return static_cast<index_t>(1);
        const std::size_t middle = nonzero.size() / 2;
        std::nth_element(nonzero.begin(), nonzero.begin() + middle, nonzero.end());
        const index_t median = nonzero[middle];
        return static_cast<index_t>(std::max<double>(
            1.0, std::ceil(std::sqrt(static_cast<double>(median)))));
    };

    left_min_candidates = adaptive_minimum(left_valid_results);
    right_min_candidates = adaptive_minimum(right_valid_results);
    const bool user_minimum = args.pdr_min_candidates > 0;
    if (user_minimum) {
        left_min_candidates = static_cast<index_t>(args.pdr_min_candidates);
        right_min_candidates = static_cast<index_t>(args.pdr_min_candidates);
    }

    auto print_stats = [&](const char* label, const std::vector<index_t>& counts) {
        std::vector<index_t> sorted = counts;
        std::sort(sorted.begin(), sorted.end());
        const index_t minimum = sorted.front();
        const index_t median = sorted[sorted.size() / 2];
        const index_t maximum = sorted.back();
        const double average = static_cast<double>(
            std::accumulate(sorted.begin(), sorted.end(), 0ULL)) /
            static_cast<double>(sorted.size());
        const std::size_t zeros = std::count(sorted.begin(), sorted.end(), 0);
        std::cout << "  " << label << " min/median/avg/max: "
                  << minimum << " / " << median << " / " << average
                  << " / " << maximum << "  | zero: " << zeros << "\n";
    };

    std::cout << "Primer3 accepted " << retval->fwd.num_elem
              << " left and " << retval->rev.num_elem
              << " right oligos across " << num_pdrs << " possible PDRs\n";
    print_stats("LEFT ", left_valid_results);
    print_stats("RIGHT", right_valid_results);
    std::cout << "  " << (user_minimum ? "User" : "Adaptive")
              << " feasibility minimum: left >= "
              << left_min_candidates << ", right >= "
              << right_min_candidates << " candidates\n";
    std::cout << "\n";

    destroy_p3retval(retval);
    destroy_seq_args(sa);
    p3_destroy_global_settings(pa);
}

bool RiskOptimizer::left_valid(index_t p) const {
    return p < left_valid_results.size() &&
           left_valid_results[p] >= left_min_candidates;
}

bool RiskOptimizer::right_valid(index_t p) const {
    return p < right_valid_results.size() &&
           right_valid_results[p] >= right_min_candidates;
}

RiskOptimizer::~RiskOptimizer() {
    delete [] risk;
    delete [] prefix_sum;
    if (max == min) {
        delete [] memo;
        delete [] prev;
    }
    delete [] solutions;
    delete [] solution_vector;
}

std::vector<index_t> RiskOptimizer::random_PDR() {
    std::vector<index_t> PDR;
    while (true) {
        index_t k = PDR.size();
        index_t cmin, cmax;
        if (k % 2 == 0) {
            if (k == 0) {
                cmin = 0;
                cmax = 3 * len - 1;
            }
            else if (k == 2) {
                if (max <= PDR[k - 1] + 3 * len)
                    cmin = std::max(PDR[k - 2] + len, PDR[k - 1] + 3 * len - max);
                else 
                    cmin = PDR[k - 2] + len;
                cmax = PDR[k - 1] - 1;
            }
            else {
                if (max <= PDR[k - 1] + 3 * len)
                    cmin = std::max(PDR[k - 3] + len, PDR[k - 1] + 3 * len - max);
                else 
                    cmin = PDR[k - 3] + len;
                cmax = PDR[k - 1] - 1;
            }
        }
        else {
            if (k == 1) {
                cmin = PDR[k - 1] + min - len;
                cmax = PDR[k - 1] + max - 1;
            }
            else {
                cmin = std::max(PDR[k - 1] - len + min, PDR[k - 2] + 2 * len);
                cmax = PDR[k - 1] + max - 1;
            }
            if (cmax >= size) {
                PDR.pop_back();
                break;
            }
        }
        PDR.push_back(greedy_random_between(cmin, cmax - len + 2));
    }
    return PDR;
}

std::vector<index_t> RiskOptimizer::random_search(std::size_t limit) {
    risk_t min = -1;
    std::vector<index_t> min_PDR;
    std::cout << std::setw(8) << "i" << std::setw(12) << "loss" << std::endl;
    for (std::size_t i = 0; i < limit; i ++) {
        std::vector<index_t> PDR = random_PDR();
        risk_t temp = score(PDR, ALPHA);
        if (min < 0 || temp < min) {
            min = temp;
            min_PDR = PDR;
            std::cout << std::setw(8) << i << std::setw(12) << min << std::endl;
        }
    }
    return min_PDR;
}

index_t RiskOptimizer::greedy_random_between(index_t cmin, index_t cmax) {
    std::vector<std::pair<index_t, risk_t>> all_sum;
    for (index_t i = cmin; i < cmax; i ++) {
        risk_t sum = 0;
        for (index_t j = i; j < i + len; j ++) 
            sum += risk[j];
        all_sum.push_back(std::make_pair(i, sum * sum));
    }
    std::sort(all_sum.begin(), all_sum.end(), [](const std::pair<index_t, risk_t> &a, const std::pair<index_t, risk_t> &b) {
        return a.second < b.second;
    });
    index_t bound = (index_t) ((cmax - cmin) * 0.3);
    if (bound == 0) bound = 1;
    return all_sum[random_between(0, bound)].first;
}

void RiskOptimizer::validate_PDR(std::vector<index_t> PDR) {
    for (index_t i = 0; i < PDR.size(); i ++) {
        const bool has_candidates =
            (i % 2 == 0) ? left_valid(PDR[i]) : right_valid(PDR[i]);
        if (!has_candidates) {
            std::cout << i << " " << PDR[i] << "\n";
        }
        assert(has_candidates);
        if (i % 2 == 0) {
            if (i == 0) {
                continue;
            }
            else if (i == 2) {
                assert(PDR[i] - PDR[i - 2] >= len);
            }
            else {
                assert(PDR[i] - PDR[i - 3] >= len);
                
            }
            index_t l = PDR[i + 1] - PDR[i];
            assert(l + len <= max);
            assert(l + len >= min);
        }
        else {
            if (i == PDR.size() - 1) {
                assert(PDR[i] - PDR[i - 2] >= len);
            }
            else {
                assert(PDR[i] - PDR[i + 1] >= len);
            }
        }

    }
    std::cout << "PDRs: " << PDR.size() << std::endl;
    std::cout << "uncovered: " << PDR[0] << " (front), " << size - PDR[PDR.size() - 1] - len << " (rear)" << std::endl;
    assert(PDR[0] <= len * 2);
    // assert(size - PDR[PDR.size() - 1] <= 3 * len);
    assert(size - PDR[PDR.size() - 1] >= len);
}

risk_t RiskOptimizer::score(std::vector<index_t> PDR, risk_t alpha) {
    std::vector<risk_t> all_sum;
    for (index_t i = 0; i < PDR.size(); i ++) {
        risk_t sum = 0;
        for (index_t j = PDR[i]; j < PDR[i] + len; j ++) 
            sum += risk[j];
        all_sum.push_back(sum * sum);
    }
    std::sort(all_sum.begin(), all_sum.end());
    risk_t total = 0;
    for (index_t i = 0; i < std::ceil(alpha * PDR.size()); i ++) 
        total += all_sum[PDR.size() - i - 1];
    return total;
}

risk_t RiskOptimizer::opt(index_t f, index_t r, index_t f_, index_t r_, risk_t u, risk_t alpha) {
    if (r < r_ + len || r_ < f + len || f < f_ + len) return INF;
    key_t k = to_key(f, r, f_, r_);
    if (memo_umap.find(k) == memo_umap.end()) {
        risk_t min_risk;
        key_t min_key;
        if (f_ <= len * 2) {
            min_risk = cost(f, u, alpha) + cost(r, u, alpha) + cost(f_, u, alpha) + cost(r_, u, alpha);
            min_key = 0;
        }
        else {
            min_risk = INF;
            for (index_t r__ = f_ + len; f >= r__ + len; r__ ++) {
                for (index_t f__ = (r__ >= max ? r__ - max : 0); r__  >= f__ + min; f__ ++) {
                    risk_t temp = opt(f_, r_, f__, r__, u, alpha);
                    if (temp < min_risk) {
                        min_risk = temp;
                        min_key = to_key(f_, r_, f__, r__);
                    }
                }
            }
            min_risk += cost(f, u, alpha) + cost(r, u, alpha);
        }
        memo_umap[k] = min_risk;
        prev_umap[k] = min_key;
    }
    return memo_umap[k];
}

risk_t RiskOptimizer::top_k_opt(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha) {
    min_PDR.clear();
    memo_umap.clear();
    prev_umap.clear();
    risk_t min_risk = INF;
    key_t min_key = 0;
    for (index_t r = size - len * 3; r <= size - len; r ++) {
        for (index_t f = (r >= max ? r - max : 0); r >= f + min; f ++) {
            for (index_t r_ = f + len; r >= r_ + len; r_ ++) {
                for (index_t f_ = (r_ >= max ? r_ - max : 0); r_ >= f_ + min; f_ ++) {
                    risk_t temp = opt(f, r, f_, r_, u, alpha);
                    if (temp < min_risk) {
                        min_risk = temp;
                        min_key = to_key(f, r, f_, r_);
                    }
                }
            }
        }
    }
    while (min_key != 0) {
        index_t a, b, c, d;
        to_index(min_key, a, b, c, d);
        min_PDR.push_back(b);
        min_PDR.push_back(a);
        min_key = prev_umap[min_key];
        if (min_key == 0) {
            min_PDR.push_back(d);
            min_PDR.push_back(c);
        }
    }
    std::reverse(min_PDR.begin(), min_PDR.end());
    return min_risk;
}

risk_t RiskOptimizer::top_k_opt_fast(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha) {
    min_PDR.clear();
    for (index_t i = 0; i < size; i ++) 
        solution_vector[i].clear();
    risk_t all_min_risk = INF;
    key_t all_min_key = 0;
    for (index_t r = max; r <= size - len; r ++) {
        index_t f = r - max;
        for (index_t r_ = f + len; r_ <= r - len; r_ ++) {
            if (r_ < max) continue;
            index_t f_ = r_ - max;
            if (r < r_ + len || r_ < f + len || f < f_ + len) continue;
            if (!right_valid(r) || !right_valid(r_) ||
                !left_valid(f) || !left_valid(f_)) continue;
            key_t k = to_key_2(r, r_);
            risk_t min_risk;
            key_t min_key;
            if (f_ <= fumax) {
                min_risk = cost(f, u, alpha) + cost(r, u, alpha) + cost(f_, u, alpha) + cost(r_, u, alpha);
                min_key = 0;
            }
            else {
                min_risk = INF;
                index_t r__min = f_ + len, r__max = f - len;
                if (r__min < max) r__min = max;
                if (r__min < f_ + len) r__min = f_ + len;
                if (r__max > r_ - len) r__max = r_ - len;
                if (r__min <= r__max) {
                    std::pair<index_t, risk_t> temp = solutions[r_]->range_query(r__min, r__max);
                    min_risk = temp.second;
                    min_key = to_key_2(r_, temp.first);
                }
                /*
                for (auto elem : solution_vector[r_]) {
                    index_t r__ = elem.first;
                    if (r__min <= r__ && r__ <= r__max) {
                        risk_t temp_risk = elem.second;
                        key_t temp_key = to_key_2(r_, r__);
                        if (temp_risk < min_risk) {
                            min_risk = temp_risk;
                            min_key = temp_key;
                        }
                    }
                }
                */
                min_risk += cost(f, u, alpha) + cost(r, u, alpha);
            }
            solution_vector[r].push_back(std::make_pair(r_, min_risk));
            prev[k] = min_key;
            if (r >= size - rumax) {
                if (min_risk < all_min_risk) {
                    all_min_risk = min_risk;
                    all_min_key = k;
                }
            }
        }
        solutions[r] = new BST(solution_vector[r]);
    }
    risk_t min_risk = all_min_risk;
    key_t min_key = all_min_key;
    while (min_key != 0) {
        index_t a, b, c, d;
        to_index_2(min_key, b, d);
        a = b - max;
        c = d - max;
        min_PDR.push_back(b);
        min_PDR.push_back(a);
        min_key = prev[min_key];
        if (min_key == 0) {
            min_PDR.push_back(d);
            min_PDR.push_back(c);
        }
    }
    std::reverse(min_PDR.begin(), min_PDR.end());
    index_t i = min_PDR.size();
    for (index_t r = max; r <= size - len; r ++) 
        delete solutions[r];
    return min_risk;
}

risk_t RiskOptimizer::opt_m(index_t f, index_t r, risk_t u, risk_t alpha) {
    key_t k = to_key_2(r, f);
    if (memo_umap.find(k) == memo_umap.end()) {
        risk_t min_risk = INF;
        key_t min_key;
        if (f <= len * 2) {
            min_risk = min_key = 0;
        }
        else {
            index_t r_min = std::max(max, f + len), r_max = (f + r + len) / 2 - len;
            for (index_t r_ = r_min; r_ <= r_max; r_ ++) {
                index_t f_min = (r_ >= max ? r_ - max : 0);
                index_t f_max = std::min((r_ >= min ? r_ - min : 0), (2 * f >= len + r_? 2 * f - len - r_ : 0));
                for (index_t f_ = f_min; f_ <= f_max; f_ ++) {
                    key_t k_ = to_key_2(r_, f_);
                    risk_t temp = opt_m(f_, r_, u, alpha);
                    if (temp < min_risk) {
                        min_risk = temp;
                        min_key = k_;
                    }
                }
            }
        }
        min_risk += cost(f, u, alpha) + cost(r, u, alpha);
        memo_umap[k] = min_risk;
        prev_umap[k] = min_key;
    }
    return memo_umap[k];
}

risk_t RiskOptimizer::top_k_opt_m(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha) {
    min_PDR.clear();
    memo_umap.clear();
    prev_umap.clear();
    risk_t min_risk = INF;
    key_t min_key = 0;
    for (index_t r = size - len * 3; r <= size - len; r ++) {
        for (index_t f = (r >= max ? r - max : 0); r >= f + min; f ++) {
            risk_t temp = opt_m(f, r, u, alpha);
            if (temp < min_risk) {
                min_risk = temp;
                min_key = to_key_2(r, f);
            }
        }
    }
    while (min_key != 0) {
        index_t a, b;
        to_index_2(min_key, b, a);
        min_PDR.push_back(b);
        min_PDR.push_back(a);
        min_key = prev_umap[min_key];
    }
    std::reverse(min_PDR.begin(), min_PDR.end());
    return min_risk;
}

risk_t RiskOptimizer::top_k_opt_mi(risk_t u, std::vector<index_t> &min_PDR, risk_t alpha) {
    min_PDR.clear();
    risk_t all_min_risk = INF;
    key_t all_min_key = 0;
    for (index_t r = max; r <= size - len; r ++) {
        for (index_t f = (r >= max ? r - max : 0); r >= f + min; f ++) {
            risk_t min_risk = INF;
            key_t min_key, k = to_key_2(r, f);
            if (f <= len * 2) {
                min_risk = min_key = 0;
            }
            else {
                index_t r_min = std::max(max, f + len), r_max = (f + r + len) / 2 - len;
                for (index_t r_ = r_min; r_ <= r_max; r_ ++) {
                    index_t f_min = (r_ >= max ? r_ - max : 0);
                    index_t f_max = std::min((r_ >= min ? r_ - min : 0), (2 * f >= len + r_? 2 * f - len - r_ : 0));
                    for (index_t f_ = f_min; f_ <= f_max; f_ ++) {
                        key_t k_ = to_key_2(r_, f_);
                        risk_t temp = memo[k_];
                        if (temp < min_risk) {
                            min_risk = temp;
                            min_key = k_;
                        }
                    }
                }
            }
            min_risk += cost(f, u, alpha) + cost(r, u, alpha);
            memo[k] = min_risk;
            prev[k] = min_key;
            if (r >= size - len * 3) {
                if (min_risk < all_min_risk) {
                    all_min_risk = min_risk;
                    all_min_key = k;
                }
            }
        }
    }
    risk_t min_risk = all_min_risk;
    key_t min_key = all_min_key;
    while (min_key != 0) {
        index_t a, b;
        to_index_2(min_key, b, a);
        min_PDR.push_back(b);
        min_PDR.push_back(a);
        min_key = prev[min_key];
    }
    std::reverse(min_PDR.begin(), min_PDR.end());
    return min_risk;
}

risk_t RiskOptimizer::top_k_opt_mi_fast(risk_t u,
                                    std::vector<index_t>& min_PDR,
                                    risk_t alpha) {
    min_PDR.clear();

    const index_t r_first = max;
    const index_t r_last = size >= len ? size - len : 0;
    if (r_first > r_last) return INF;

    if (!midpoint_rmq_) {
        midpoint_rmq_ = std::make_unique<MidpointRMQ2D<risk_t, key_t>>(
            r_first, r_last, min, max);
    }
    midpoint_rmq_->reset();

    risk_t all_min_risk = INF;
    key_t  all_min_key  = 0;

    const index_t terminal_r = size > rumax ? size - rumax : 0;
    for (index_t r = r_first; r <= r_last; ++r) {
        const index_t f_first = r - max;
        const index_t f_last = r - min;

        for (index_t f = f_first; f <= f_last; ++f) {
            if (!right_valid(r) || !left_valid(f)) continue;

            risk_t min_risk = INF;
            key_t  min_key  = 0;
            key_t  k        = to_key_2(r, f);

            if (f <= fumax) {
                min_risk = 0;
            } else {
                const index_t pred_r_min = std::max(max, f + len);
                const index_t half_sum = (f + r + len) / 2;
                if (half_sum >= len && 2 * f >= len) {
                    const index_t pred_r_max = half_sum - len;
                    const index_t midpoint_sum_max = 2 * f - len;
                    const auto best = midpoint_rmq_->query(
                        pred_r_min, pred_r_max, midpoint_sum_max);
                    if (best.valid) {
                        min_risk = best.value;
                        min_key = best.key;
                    }
                }
            }

            if (min_risk >= INF) continue;
            min_risk += cost(f, u, alpha) + cost(r, u, alpha);
            memo[k]   = min_risk;
            prev[k]   = min_key;

            midpoint_rmq_->insert(r, r + f, min_risk, k);

            if (r >= terminal_r) {
                if (min_risk < all_min_risk) {
                    all_min_risk = min_risk;
                    all_min_key  = k;
                }
            }
        }
    }

    // traceback
    key_t min_key = all_min_key;
    while (min_key != 0) {
        index_t a, b;
        to_index_2(min_key, b, a);
        min_PDR.push_back(b);
        min_PDR.push_back(a);
        min_key = prev[min_key];
    }
    std::reverse(min_PDR.begin(), min_PDR.end());
    return all_min_risk;
}

std::vector<index_t> RiskOptimizer::search(risk_t lower, 
                                           risk_t upper, 
                                           risk_t eps) {
    const int TW = 8 + 12 + 12 + 12; // = 44
    std::cout << "Search for optimal PDR...\n"
              << std::string(TW, '-') << "\n"
              << std::right
              << std::setw(8)  << "iter"
              << std::setw(12) << "u'"
              << std::setw(12) << "loss'"
              << std::setw(12) << "loss"
              << "\n"
              << std::string(TW, '-') << "\n";

    max -= len; min -= len;
    midpoint_rmq_.reset();
    std::vector<index_t> min_PDR, final_PDR, temp_PDR;
    risk_t alpha = ALPHA, l = lower, r = upper;

    // A global count statistic does not guarantee that a complete alternating
    // left/right PDR chain is feasible at that minimum. Test graph feasibility
    // once before the numerical search. Adaptive minima may be relaxed, but an
    // explicit user minimum remains a strict constraint.
    auto find_feasible_chain = [&]() {
        min_PDR.clear();
        const risk_t probe_u = (lower + upper) / 2;
        if (min == max)
            top_k_opt_fast(probe_u, min_PDR, alpha);
        else
            top_k_opt_mi_fast(probe_u, min_PDR, alpha);
        return !min_PDR.empty();
    };

    bool feasible = find_feasible_chain();
    if (!feasible && args.pdr_min_candidates == 0) {
        const index_t initial_left = left_min_candidates;
        const index_t initial_right = right_min_candidates;
        while (!feasible &&
               (left_min_candidates > 1 || right_min_candidates > 1)) {
            if (left_min_candidates > 1) --left_min_candidates;
            if (right_min_candidates > 1) --right_min_candidates;
            feasible = find_feasible_chain();
        }
        if (feasible) {
            std::cout
                << "Adaptive feasibility minimum relaxed for a complete PDR "
                   "chain: left "
                << initial_left << " -> " << left_min_candidates
                << ", right " << initial_right << " -> "
                << right_min_candidates << "\n\n";
        }
    }

    if (!feasible) {
        max += len;
        min += len;
        std::string message =
            "no complete PDR chain satisfies the Primer3 candidate "
            "constraints (left >= " + std::to_string(left_min_candidates) +
            ", right >= " + std::to_string(right_min_candidates) + ")";
        if (args.pdr_min_candidates > 0) {
            message +=
                "; reduce --pdr-min-candidates or use 0 for adaptive mode";
        } else {
            message +=
                "; no complete chain exists even with one candidate per side";
        }
        throw std::runtime_error(message);
    }

    min_PDR.clear();

    for (index_t i = 1; r - l > eps; i++) {
        risk_t u = (l + r) / 2, temp;
        if (min == max)
            temp = top_k_opt_fast(u, min_PDR, alpha);
        else
            temp = top_k_opt_mi_fast(u, min_PDR, alpha);

        if (final_PDR.size() == 0 || score(min_PDR, alpha) < score(final_PDR, alpha))
            final_PDR = min_PDR;

        risk_t grad = alpha * min_PDR.size();
        for (index_t j = 0; j < min_PDR.size(); j++)
            if (u < cost(min_PDR[j], 0, 0)) grad--;
        (grad < 0 ? l : r) = u;

        std::cout << std::right
                  << std::setw(8)  << i
                  << std::setw(12) << u
                  << std::setw(12) << temp
                  << std::setw(12) << score(min_PDR, alpha)
                  << "\n";
    }

    std::cout << std::string(TW, '-') << "\n";
    max += len; min += len;
    std::cout << "Search completed!\n\n";
    return final_PDR;
}

