#pragma once
#include "utility.hpp"

struct PDRLayoutWindow {
    std::int64_t left_start, left_end, right_start, right_end;
};

// Empty result means valid. Checks geometry, not primer feasibility or scores.
std::string pdr_layout_error(const std::vector<PDRLayoutWindow>& windows,
                            std::size_t target_length, std::size_t front_limit,
                            std::size_t rear_limit);

// Expansion geometry: alternating pools, disjoint same-pool product envelopes,
// and overlapping adjacent product envelopes for every possible primer choice.
// Unlike initial midpoint placement, no midpoint anchors are imposed.
std::string expanded_pdr_layout_error(const std::vector<PDRLayoutWindow>& windows,
                                     std::size_t target_length, std::size_t front_limit,
                                     std::size_t rear_limit);

// Expand only affected primer windows, preferring symmetric growth. If a side
// is blocked, the other may grow, subject to the complete expanded layout.
void expand_pdr_window(std::vector<PDRLayoutWindow>& windows, std::size_t index,
                       bool left_affected, bool right_affected, std::size_t bases,
                       std::size_t target_length, std::size_t front_limit,
                       std::size_t rear_limit);
