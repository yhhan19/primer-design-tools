#include "pdr_layout.hpp"

std::string pdr_layout_error(const std::vector<PDRLayoutWindow>& windows,
                            std::size_t target_length, std::size_t front_limit,
                            std::size_t rear_limit) {
    if (windows.empty()) return "empty PDR layout";
    for (std::size_t i = 0; i < windows.size(); ++i) {
        const auto& w = windows[i];
        if (w.left_start < 0 || w.left_end <= w.left_start ||
            w.left_end > std::int64_t(target_length) ||
            w.right_start < w.left_start || w.right_end <= w.right_start ||
            w.right_end > std::int64_t(target_length))
            return "PDR window bounds/order at amplicon " + std::to_string(i);
        if (i) {
            const auto& p = windows[i-1];
            if (p.left_start + p.right_end > 2 * w.left_start)
                return "forward PDR precedes previous midpoint at amplicon " + std::to_string(i);
            if (w.left_end > p.right_start)
                return "forward/reverse PDR interleaving at amplicon " + std::to_string(i);
            if (2 * p.right_end > w.left_start + w.right_end)
                return "previous reverse PDR exceeds current midpoint at amplicon " + std::to_string(i);
        }
    }
    if (windows.front().left_start > std::int64_t(front_limit)) return "front PDR limit";
    const auto terminal = target_length > rear_limit ? target_length - rear_limit : 0;
    if (windows.back().right_start < std::int64_t(terminal)) return "rear PDR limit";
    return "";
}

std::string expanded_pdr_layout_error(const std::vector<PDRLayoutWindow>& windows,
                                     std::size_t target_length, std::size_t front_limit,
                                     std::size_t rear_limit) {
    if (windows.empty()) return "empty PDR layout";
    for (std::size_t i = 0; i < windows.size(); ++i) {
        const auto& w = windows[i];
        if (w.left_start < 0 || w.left_end <= w.left_start ||
            w.left_end > std::int64_t(target_length) || w.right_start < w.left_start ||
            w.right_end <= w.right_start || w.right_end > std::int64_t(target_length))
            return "expanded PDR bounds/order at amplicon " + std::to_string(i);
        // The latest possible start of this product must not leave a gap
        // after the earliest possible end of the preceding product.
        if (i && w.left_end > windows[i-1].right_start)
            return "adjacent amplicon coverage at amplicon " + std::to_string(i);
        // Alternating pool assignments: i and i-2 belong to the same pool.
        // Half-open product envelopes may touch, but cannot overlap.
        if (i >= 2 && windows[i-2].right_end > w.left_start)
            return "same-pool product envelopes overlap at amplicon " + std::to_string(i);
    }
    if (windows.front().left_start > std::int64_t(front_limit)) return "front PDR limit";
    const auto terminal = target_length > rear_limit ? target_length - rear_limit : 0;
    if (windows.back().right_start < std::int64_t(terminal)) return "rear PDR limit";
    return "";
}

void expand_pdr_window(std::vector<PDRLayoutWindow>& windows, std::size_t index,
                       bool left_affected, bool right_affected, std::size_t bases,
                       std::size_t target_length, std::size_t front_limit,
                       std::size_t rear_limit) {
    if (index >= windows.size() || !expanded_pdr_layout_error(windows, target_length, front_limit, rear_limit).empty())
        throw std::invalid_argument("Expansion requires a valid alternating-pool PDR layout");
    if (!left_affected && !right_affected) return;
    const auto original = windows[index];
    const auto limit = std::min(bases, target_length);
    auto try_change = [&](const PDRLayoutWindow& proposed) {
        const auto old = windows[index];
        windows[index] = proposed;
        if (!expanded_pdr_layout_error(windows, target_length, front_limit, rear_limit).empty()) {
            windows[index] = old;
            return false;
        }
        return true;
    };
    for (std::size_t step = 0; step < limit; ++step) {
        bool changed = false;
        for (bool left : {true, false}) {
            if (!(left ? left_affected : right_affected)) continue;
            auto w = windows[index];
            auto& start = left ? w.left_start : w.right_start;
            auto& end = left ? w.left_end : w.right_end;
            --start; ++end;
            if (try_change(w)) { changed = true; continue; }
            w = windows[index];
            --(left ? w.left_start : w.right_start);
            changed |= try_change(w);
            w = windows[index];
            ++(left ? w.left_end : w.right_end);
            changed |= try_change(w);
        }
        if (!changed) break;
    }
    const auto& w = windows[index];
    if (w.left_start > original.left_start || w.left_end < original.left_end ||
        w.right_start > original.right_start || w.right_end < original.right_end)
        throw std::logic_error("Expansion shrank an original PDR");
}
