#pragma once

#include <cstring>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <chrono>
#include <algorithm>
#include <condition_variable>
#include <cstddef>
#include <deque>
#include <mutex>
#include <thread>
#include <stdexcept>
#include <cmath>
#include <ctime>
#include <cctype>
#include <iomanip>
#include <iterator>
#include <limits>
#include <numeric>
#include <atomic>
#include <set>
#include <map>
#include <future>
#include <filesystem>
#include <memory>
#include <optional>

namespace fs = std::filesystem;

extern "C" {
#include "thal.h"
#include "libprimer3.h"
#include "p3_seq_lib.h"
#include "oligotm.h"
}

#define key_t    uint64_t
#define index_t  uint32_t
#define risk_t   double
#define weight_t double

#define ALPHA       0.1
#define INF         1e18

#define ALPHABET_SIZE 128
#define KELVIN_ZERO   -273.15

#ifndef PRIMER3_PATH
#define PRIMER3_PATH "./primer3"
#endif

struct DimerPartner {
    std::string primer_name;
    double dg;
};

struct PrimerDimerRisk {
    std::optional<double> cutoff;
    std::vector<DimerPartner> partners;
};

struct Oligo {
    std::string seq;
    int start;
    int length;
    double tm;
    double gc;

    PrimerDimerRisk dimer_risk;
};

struct PairMsaMetrics {
    double joint_coverage = 0;
    double joint_perfect_match = 0;
    double joint_3prime_match = 0;
    double weighted_mismatch = 1;
    double unknown_fraction = 1;
    std::size_t sequence_count = 0;
};

struct PrimerResult {
    PairMsaMetrics msa;
    Oligo       left, right;
    std::size_t product_size = 0;
    index_t     pdr_left     = 0;
    index_t     pdr_right    = 0;
    std::size_t segment_id   = 0;
    std::size_t pool         = 0;
};

struct PrimerOutput {
    bool pairs_required = false;  // When set, only combinations in pairs may be selected.
    std::vector<PrimerResult> pairs;
    std::vector<Oligo>        left_oligos;
    std::vector<Oligo>        right_oligos;
    index_t     pdr_left     = 0;
    index_t     pdr_right    = 0;
    std::size_t segment_id  = 0;
};

struct Result {
    std::string label;
    std::string data;
    std::string ref_label;
    std::string ref_seq;
    std::size_t ref_index;
    std::size_t flag;
    double      dg;
};

struct Args {
    std::string log_file;
    std::string console_mode = "auto";
    // std::string mode;

    std::string input_file;
    std::string output_file;
    std::string ref_file;
    std::string rate_file;

    std::size_t len_amp     = 420;
    std::size_t len_amp_min = 252;
    std::size_t len_PDR     = 40;
    std::size_t pdr_min_candidates = 0; // 0 = adaptive
    std::string pdr_search = "midpoint";  // fixed or midpoint
    std::size_t len_ru      = 80;
    std::size_t len_fu      = 80;
    double      u_max       = 10000;
    double      u_min       = 0.1;

    std::size_t iter        = 1000;
    std::size_t dimer_restarts = 400;
    std::size_t seed        = 42;

    std::string off_target_mode = "paired";
    std::size_t off_target_amp_min = 1;
    std::size_t off_target_amp_max = 1000;
    std::size_t kmer_len    = 8;
    std::size_t threshold   = 4;
    std::size_t nthreads    = 8;
    std::size_t block_size  = 1 << 27;
    std::size_t chunk_size  = (std::size_t) 1e8;

    std::size_t p3_opt_size = 20;
    std::size_t p3_min_size = 18;
    std::size_t p3_max_size = 25;
    double      p3_opt_tm   = 60.0;
    double      p3_min_tm   = 57.0;
    double      p3_max_tm   = 63.0;
    double      p3_min_gc   = 40.0;
    double      p3_max_gc   = 60.0;
    std::size_t num_return  = 5;
    double pair_keep_percent = 20.0;
    std::size_t pair_min_return = 5;
    std::size_t pair_max_mismatches = 2;
    std::size_t pair_max_3prime_mismatches = 0;
    std::size_t pair_3prime_window = 5;

    double      mv          = 50.0;
    double      dv          = 1.5;
    double      dntp        = 0.6;
    double      dna_conc    = 200.0;
    
    double      dg_thres    = -15000;
    double      temp        = 37.0;
    // Dimer window defaults to 58-62 C, independently of non-dimer temp.
    static constexpr double default_dimer_dg_cutoff = -3000.0;  // cal/mol
    std::string dimer_methods = "min-affected-primers";
    bool expand_affected_pdrs = true;
    std::size_t pdr_expansion_bases = 20;
    std::size_t dimer_rounds = 1;
    std::size_t dimer_landscape_points = 0;
    bool dimer_landscape_adaptive = false;
    bool landscape_only = false;
    std::size_t dimer_landscape_restarts = 200;
    std::size_t dimer_landscape_iterations = 1000;
    double dimer_landscape_min = -5000;
    double dimer_landscape_max = -2000;
    std::string dimer_start_prefix;
    std::string dimer_cache = "auto";
    std::string dimer_select;  // Empty selects the first listed method.
    std::optional<double> dimer_dg_cutoff;  // cal/mol
    double dimer_bad_percent = 1.0;  // Percentage of selected pool edges allowed below the rank threshold.
    std::optional<double> dimer_temp_min = 58.0;
    std::optional<double> dimer_temp_max = 62.0;

    bool        help        = false;
    bool help_all = false;
};

struct SegmentContext {
    std::string               name;
    std::vector<std::string>  labels, sequences;
    std::string               tmpl;
    std::vector<index_t>      pdr_regions;
    std::vector<std::size_t> pdr_region_lengths;
    std::vector<PrimerOutput> candidate_primers, filtered_primers;
    std::vector<Result>       off_target_hits;
};

struct PipelineContext {
    Args                         args;
    std::vector<SegmentContext>  segments;

    // flat members — stages use these directly, unchanged
    std::string               tmpl;
    std::vector<std::string>  labels, sequences;
    std::vector<index_t>      pdr_regions;
    std::vector<PrimerOutput> candidate_primers, filtered_primers;
    std::vector<Result>       off_target_hits;

    // global
    std::vector<index_t>      dimer_solution;
    std::vector<std::size_t> pdr_region_lengths; // Optional per-window generation lengths.
    std::vector<PrimerResult> solution_primers;

    void load_segment(std::size_t i) {
        current_seg_      = i;
        auto& seg         = segments[i];
        tmpl              = seg.tmpl;
        labels            = seg.labels;
        sequences         = seg.sequences;
        pdr_regions       = seg.pdr_regions;
        pdr_region_lengths = seg.pdr_region_lengths;
        candidate_primers = seg.candidate_primers;
        filtered_primers  = seg.filtered_primers;
        off_target_hits   = seg.off_target_hits;
    }

    void save_segment() {
        auto& seg             = segments[current_seg_];
        seg.tmpl              = tmpl;
        seg.labels            = labels;
        seg.sequences         = sequences;
        seg.pdr_regions       = pdr_regions;
        seg.pdr_region_lengths = pdr_region_lengths;
        seg.candidate_primers = candidate_primers;
        seg.filtered_primers  = filtered_primers;
        seg.off_target_hits   = off_target_hits;
    }

    std::size_t current_seg_id() const { return current_seg_; }

private:
    std::size_t current_seg_ = 0;
};

std::size_t read_fasta(const std::string& file_name,
                       std::vector<std::string>& labels,
                       std::vector<std::string>& data);

void write_fasta(const std::string& filename,
                 const std::vector<std::string>& headers,
                 const std::vector<std::string>& sequences,
                 int line_width);

void char_freq(const std::string& seq);

std::size_t read_rate(std::string filename,
                      std::vector<risk_t>& input);

std::size_t read_rate_2(std::string filename,
                        std::vector<risk_t>& input,
                        std::string& ref);

std::size_t read_rate_ref(std::string filename,
                          std::vector<risk_t>& input,
                          std::string& ref);

void write_results(const std::string& output_file,
                   const std::vector<Result>& results,
                   const std::vector<std::string>& labels);

void test_primer3();

key_t to_key  (index_t f, index_t r, index_t f_, index_t r_);
void  to_index(key_t k, index_t& f, index_t& r, index_t& f_, index_t& r_);

key_t to_key_2  (index_t r, index_t r_);
void  to_index_2(key_t k, index_t& r, index_t& r_);

bool        is_ws  (char c);
char        comp   (char c);
std::string revcomp(const std::string& s);
bool        is_gc  (char c);

std::string display_PDR(std::vector<index_t> PDR);

std::unordered_map<std::string, std::vector<Result>>
group_by_label(const std::vector<Result>& results);

std::vector<risk_t> random_risk   (std::size_t size);
std::size_t         random_between(std::size_t min, std::size_t max);

std::string msa_consensus(const std::vector<std::string>& msa);

void display_primer_output(const PrimerOutput& out);

void convert(const std::vector<PrimerOutput>& primer_output,
             std::vector<std::string>& labels,
             std::vector<std::string>& sequences);


Oligo extract_oligo(const primer_rec& p,
                    const char* sequence,
                    bool is_right);

PrimerOutput extract_all(const p3retval* retval,
                         const seq_args*  sa);

void print_solution(const std::vector<PrimerResult>& solution);

void save_solution(std::vector<SegmentContext>& seg, std::vector<PrimerResult>& solution, const std::string& filename);
void annotate_solution_dimers(const std::vector<SegmentContext>& segments,
                              std::vector<PrimerResult>& solution,
                              double cutoff, double min_c, double max_c,
                              const std::string& output_prefix);

void save_coverage_report(const std::vector<SegmentContext>& segments,
                          const std::vector<PrimerResult>& solution,
                          const std::string& output_prefix);

bool same_oligo(const Oligo& a, const Oligo& b);
struct PairingStats {
    std::size_t combinations = 0;
    std::size_t accepted = 0;
};
PairingStats build_candidate_pairs(PrimerOutput& output, const std::string& target, const Args& args);
void retain_surviving_pairs(PrimerOutput& output);
void score_candidate_pairs_msa(PrimerOutput& output,
    const std::vector<std::string>& sequences, const Args& args);
std::size_t retain_candidate_pairs_msa(PrimerOutput& output, const Args& args);

struct RemovalSets {
    std::set<int> remove_pair_F, remove_pair_R, remove_left, remove_right;
    std::size_t num_return;
    bool all_removed(const PrimerOutput& original) const {
        size_t kept_left  = original.left_oligos.size()  - remove_left.size();
        size_t kept_right = original.right_oligos.size() - remove_right.size();
        if (original.pairs_required) {
            bool any_pair = false;
            for (const auto& pair : original.pairs) {
                bool left = false, right = false;
                for (std::size_t i = 0; i < original.left_oligos.size(); ++i)
                    if (!remove_left.count(i) && same_oligo(pair.left, original.left_oligos[i])) left = true;
                for (std::size_t i = 0; i < original.right_oligos.size(); ++i)
                    if (!remove_right.count(i) && same_oligo(pair.right, original.right_oligos[i])) right = true;
                any_pair = any_pair || (left && right);
            }
            if (!any_pair) return true;
        }
        return kept_left < num_return || kept_right < num_return;
    }
};

RemovalSets computeRemoval(const std::vector<const Result*>& output_results,
                           double threshold, 
                           std::size_t num_return) ;

std::vector<PrimerOutput> filterByDG_relax(
    const std::vector<Result>& results,
    double dg_threshold,
    const std::vector<PrimerOutput>& original_primers,
    std::size_t num_return);
