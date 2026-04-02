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
#include <numeric>
extern "C" {
#include "thal.h"
#include "libprimer3.h"
#include "p3_seq_lib.h"
#include "oligotm.h"
}

#define risk_t double
#define ALPHA 0.1
#define key_t uint64_t
#define index_t uint32_t
#define INF 1e18
#define ITER_LIMIT 500
#define KEY_LIMIT (1 << 26)
#define weight_t double

std::vector<risk_t> random_risk(std::size_t size);
std::size_t random_between(std::size_t min, std::size_t max);
key_t to_key(index_t f, index_t r, index_t f_, index_t r_);
void to_index(key_t k, index_t &f, index_t &r, index_t &f_, index_t &r_);
std::size_t read_rate(std::string filename, std::vector<risk_t> &input);
key_t to_key_2(index_t r, index_t r_);
void to_index_2(key_t k, index_t &r, index_t &r_);
std::size_t read_rate_2(std::string filename, 
                        std::vector<risk_t> &input, 
                        std::string &ref);
std::string display_PDR(std::vector<index_t> PDR);
std::size_t read_rate_ref(std::string filename, 
                          std::vector<risk_t> &input, 
                          std::string &ref);


#define ALPHABET_SIZE 128
#define KELVIN_ZERO -273.15
#ifndef PRIMER3_PATH
#define PRIMER3_PATH "./primer3"
#endif

std::size_t read_fasta(const std::string &file_name, 
                       std::vector<std::string> &labels, 
                       std::vector<std::string> &data);

bool is_ws(char c);
char comp(char c);
std::string revcomp(const std::string& s);

struct Result {
    std::string label;
    std::string data;
    std::string ref_label;
    std::string ref_seq;
    std::size_t ref_index;
    std::size_t flag;
    double dg;
};

std::unordered_map<std::string, std::vector<Result>>
group_by_label(const std::vector<Result>& results);

void write_results(const std::string& output_file,
                   const std::vector<Result>& results,
                   const std::vector<std::string>& labels);

struct Args {
    std::string mode;

    std::string input_file;
    std::string output_file;
    std::string ref_file;

    std::size_t seed = 42;
    std::size_t len_amp = 420;
    std::size_t len_amp_min = 252;
    std::size_t len_PDR = 40;
    double u_max = 10000;
    double u_min = 0.1;

    std::size_t iter = 1000;

    std::size_t kmer_len = 8;
    std::size_t threshold = 2;
    std::size_t nthreads = 8;
    std::size_t block_size = 1 << 27;
    std::size_t chunk_size = (std::size_t) 1e8;
    
    double dg_thres = 1e8;
    double mv = 50.0;
    double dv = 1.5;
    double dntp = 0.6;
    double dna_conc = 200.0;
    double temp = 37.0;

    bool help = false;
};

struct PipelineContext {
    Args args;

    std::vector<std::string> labels, sequences;
    std::vector<index_t> pdr_regions;
    std::vector<index_t> candidate_pairs;
    std::vector<index_t> dimer_solution;
    std::vector<index_t> off_target_hits;
};
