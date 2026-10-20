#pragma once
#include <vector>
#include <limits>
#include <cstddef>
#include <algorithm>

template<typename Value, typename Key>
class BandRMQ2D {
public:
    struct Entry {
        Value value = std::numeric_limits<Value>::max();
        Key   key   = Key{};
        bool  valid = false;
    };

    // r_size:   total range of r values
    // band:     max - min (width of f relative to r)
    // r_offset: minimum value of r (e.g. = max parameter)
    BandRMQ2D() {}

    BandRMQ2D(std::size_t r_size, std::size_t band, std::size_t r_offset)
        : r_size_(r_size), band_(band), r_offset_(r_offset),
          log_r_(log2_floor(r_size) + 1),
          log_b_(log2_floor(band)   + 1),
          data_(r_size * band),
          sparse_(log_r_ * log_b_ * r_size * band)
    {}

    // f_rel = f - (r - band_max) maps f into [0, band)
    // caller passes actual r and f values
    void insert(std::size_t r, std::size_t f, Value value, Key key) {
        std::size_t ri = r - r_offset_;
        std::size_t fi = f_to_band(r, f);
        if (ri >= r_size_ || fi >= band_) return;

        auto& e = data_[idx(ri, fi)];
        if (!e.valid || value < e.value)
            e = { value, key, true };
    }

    void rebuild() {
        // level (0,0)
        for (std::size_t ri = 0; ri < r_size_; ++ri)
            for (std::size_t fi = 0; fi < band_; ++fi)
                sparse_[sidx(0, 0, ri, fi)] = data_[idx(ri, fi)];

        // extend along r
        for (std::size_t lr = 1; lr < log_r_; ++lr) {
            std::size_t sr = 1 << (lr - 1);
            for (std::size_t ri = 0; ri < r_size_; ++ri)
                for (std::size_t fi = 0; fi < band_; ++fi) {
                    Entry a = sparse_[sidx(lr-1, 0, ri, fi)];
                    Entry b = (ri + sr < r_size_)
                              ? sparse_[sidx(lr-1, 0, ri + sr, fi)]
                              : Entry{};
                    sparse_[sidx(lr, 0, ri, fi)] = min_entry(a, b);
                }
        }

        // extend along band
        for (std::size_t lr = 0; lr < log_r_; ++lr)
            for (std::size_t lb = 1; lb < log_b_; ++lb) {
                std::size_t sb = 1 << (lb - 1);
                for (std::size_t ri = 0; ri < r_size_; ++ri)
                    for (std::size_t fi = 0; fi < band_; ++fi) {
                        Entry a = sparse_[sidx(lr, lb-1, ri, fi)];
                        Entry b = (fi + sb < band_)
                                  ? sparse_[sidx(lr, lb-1, ri, fi + sb)]
                                  : Entry{};
                        sparse_[sidx(lr, lb, ri, fi)] = min_entry(a, b);
                    }
            }
    }

    // query: actual r and f values
    Entry query(std::size_t r_lo, std::size_t r_hi,
                std::size_t f_lo, std::size_t f_hi) const {
        if (r_lo > r_hi || f_lo > f_hi) return {};

        // convert to internal indices
        // fi is relative to band start — use conservative bounds
        std::size_t ri_lo = (r_lo >= r_offset_) ? r_lo - r_offset_ : 0;
        std::size_t ri_hi = std::min(r_hi - r_offset_, r_size_ - 1);

        // f band indices — fi = f - (r - band_) but r varies in query
        // use global f_lo/f_hi mapped conservatively
        std::size_t fi_lo = (f_lo < band_) ? 0       : f_lo - band_;
        std::size_t fi_hi = std::min(f_hi, band_ - 1);

        if (ri_lo > ri_hi || fi_lo > fi_hi) return {};

        std::size_t lr = log2_floor(ri_hi - ri_lo + 1);
        std::size_t lb = log2_floor(fi_hi - fi_lo + 1);
        std::size_t sr = 1 << lr;
        std::size_t sb = 1 << lb;

        Entry c[4] = {
            sparse_[sidx(lr, lb, ri_lo,          fi_lo         )],
            sparse_[sidx(lr, lb, ri_lo,          fi_hi - sb + 1)],
            sparse_[sidx(lr, lb, ri_hi - sr + 1, fi_lo         )],
            sparse_[sidx(lr, lb, ri_hi - sr + 1, fi_hi - sb + 1)],
        };

        Entry best;
        for (auto& e : c)
            if (e.valid && e.value < best.value) best = e;
        return best;
    }

private:
    std::size_t r_size_, band_, r_offset_;
    std::size_t log_r_,  log_b_;
    std::vector<Entry> data_;
    std::vector<Entry> sparse_;

    // map actual f to band index given r
    std::size_t f_to_band(std::size_t r, std::size_t f) const {
        // f ∈ [r - band_, r - ?] — fi = f - (r - band_)
        std::size_t f_base = (r >= band_) ? r - band_ : 0;
        return (f >= f_base) ? f - f_base : 0;
    }

    std::size_t idx(std::size_t ri, std::size_t fi) const {
        return ri * band_ + fi;
    }

    std::size_t sidx(std::size_t lr, std::size_t lb,
                     std::size_t ri, std::size_t fi) const {
        return ((lr * log_b_ + lb) * r_size_ + ri) * band_ + fi;
    }

    static std::size_t log2_floor(std::size_t n) {
        if (n == 0) return 0;
        std::size_t k = 0;
        while ((1u << (k+1)) <= n) ++k;
        return k;
    }

    static Entry min_entry(const Entry& a, const Entry& b) {
        if (!a.valid) return b;
        if (!b.valid) return a;
        return a.value <= b.value ? a : b;
    }
};

// Dynamic 2D range-minimum structure for midpoint-constrained PDR DP.
//
// A DP state (left=f, right=r) is stored at (r, r+f). The predecessor
// midpoint constraint r'+f' <= limit then becomes a prefix query in the
// second dimension, while the allowed r' interval is a standard range query
// in the first dimension. An outer segment tree over r contains compressed
// Fenwick trees over the possible r+f interval. Both update and query are
// O(log^2 n). The interval representation avoids duplicating every state at
// every outer-tree level.
template<typename Value, typename Key>
class MidpointRMQ2D {
public:
    struct Entry {
        Value value = std::numeric_limits<Value>::max();
        Key key = Key{};
        bool valid = false;
    };

    MidpointRMQ2D(std::size_t r_min,
                  std::size_t r_max,
                  std::size_t span_min,
                  std::size_t span_max)
        : r_min_(r_min),
          r_max_(r_max),
          leaf_count_(r_max >= r_min ? r_max - r_min + 1 : 0),
          tree_base_(1) {
        while (tree_base_ < leaf_count_) tree_base_ <<= 1;
        lower_.resize(tree_base_ * 2);
        upper_.resize(tree_base_ * 2);
        fenwick_.resize(tree_base_ * 2);

        if (leaf_count_ == 0) return;
        initialize_node(1, 0, tree_base_ - 1, span_min, span_max);
    }

    void reset() {
        for (auto& tree : fenwick_) {
            std::fill(tree.begin(), tree.end(), Entry{});
        }
    }

    void insert(std::size_t r,
                std::size_t midpoint_sum,
                Value value,
                Key key) {
        if (r < r_min_ || r > r_max_) return;
        const Entry candidate{value, key, true};
        for (std::size_t node = tree_base_ + (r - r_min_);
             node > 0; node >>= 1) {
            if (fenwick_[node].empty() || midpoint_sum < lower_[node] ||
                midpoint_sum > upper_[node]) continue;
            std::size_t index = midpoint_sum - lower_[node] + 1;
            for (; index < fenwick_[node].size();
                 index += lowbit(index)) {
                fenwick_[node][index] = min_entry(
                    fenwick_[node][index], candidate);
            }
        }
    }

    Entry query(std::size_t r_lo,
                std::size_t r_hi,
                std::size_t midpoint_sum_hi) const {
        Entry best;
        if (leaf_count_ == 0 || r_lo > r_hi ||
            r_hi < r_min_ || r_lo > r_max_) return best;

        r_lo = std::max(r_lo, r_min_);
        r_hi = std::min(r_hi, r_max_);
        std::size_t left = tree_base_ + (r_lo - r_min_);
        std::size_t right = tree_base_ + (r_hi - r_min_);
        while (left <= right) {
            if (left & 1U) best = min_entry(
                best, query_node(left++, midpoint_sum_hi));
            if (!(right & 1U)) best = min_entry(
                best, query_node(right--, midpoint_sum_hi));
            left >>= 1;
            right >>= 1;
        }
        return best;
    }

private:
    std::size_t r_min_ = 0;
    std::size_t r_max_ = 0;
    std::size_t leaf_count_ = 0;
    std::size_t tree_base_ = 1;
    std::vector<std::size_t> lower_;
    std::vector<std::size_t> upper_;
    std::vector<std::vector<Entry>> fenwick_;

    void initialize_node(std::size_t node,
                         std::size_t leaf_lo,
                         std::size_t leaf_hi,
                         std::size_t span_min,
                         std::size_t span_max) {
        if (leaf_lo >= leaf_count_) return;
        const std::size_t populated_hi =
            std::min(leaf_hi, leaf_count_ - 1);
        const std::size_t first_r = r_min_ + leaf_lo;
        const std::size_t last_r = r_min_ + populated_hi;
        lower_[node] = 2 * first_r - span_max;
        upper_[node] = 2 * last_r - span_min;
        fenwick_[node].assign(upper_[node] - lower_[node] + 2, Entry{});

        if (leaf_lo == leaf_hi) return;
        const std::size_t midpoint = (leaf_lo + leaf_hi) / 2;
        initialize_node(node * 2, leaf_lo, midpoint, span_min, span_max);
        initialize_node(node * 2 + 1, midpoint + 1, leaf_hi,
                        span_min, span_max);
    }

    Entry query_node(std::size_t node,
                     std::size_t midpoint_sum_hi) const {
        Entry best;
        if (fenwick_[node].empty() || midpoint_sum_hi < lower_[node])
            return best;
        const std::size_t capped = std::min(midpoint_sum_hi, upper_[node]);
        std::size_t index = capped - lower_[node] + 1;
        for (; index > 0; index -= lowbit(index)) {
            best = min_entry(best, fenwick_[node][index]);
        }
        return best;
    }

    static std::size_t lowbit(std::size_t value) {
        return value & (~value + 1);
    }

    static Entry min_entry(const Entry& a, const Entry& b) {
        if (!a.valid) return b;
        if (!b.valid) return a;
        return a.value <= b.value ? a : b;
    }
};
