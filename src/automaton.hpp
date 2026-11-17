#pragma once
#include "utility.hpp"
#include "task.hpp"
#include "thal.hpp"

struct CandidateKey {
    std::size_t seq_id;
    std::size_t start;

    bool operator==(const CandidateKey& other) const {
        return seq_id == other.seq_id && start == other.start;
    }
};

struct CandidateKeyHash {
    std::size_t operator()(const CandidateKey& x) const {
        return std::hash<std::size_t>{}(x.seq_id) ^
               (std::hash<std::size_t>{}(x.start) << 1);
    }
};

struct ActiveEntry {
    CandidateKey key;
    std::size_t expire_at;
};

class AutomatonNode {
public:
    AutomatonNode();
    ~AutomatonNode();

private:
    friend class Automaton;

    AutomatonNode *children[ALPHABET_SIZE], *prev, *parent;
    std::vector<std::size_t> index;
    std::size_t ch;
};

struct AmpliconPrimerPair {
    std::string first_label, second_label;
    bool first_is_right = false, second_is_right = true;
    std::size_t output_index = 0, candidate_pair_index = 0;
};
struct AmpliconHit {
    std::size_t pair_index = 0;
    std::string contig;
    std::size_t start = 0, end = 0;  // Global zero-based, half-open product interval.
    Result forward_hit, reverse_hit;  // Genomic order, not candidate F/R names.
};

class Automaton {
public:
    Automaton(const std::vector<std::string> &labels, 
              const std::vector<std::string> &seqs, 
              std::size_t k);
    ~Automaton();
    std::vector<Result> search(const std::string &path, 
                               std::size_t overlap, 
                               std::size_t threshold,
                               double dg_thres,
                               std::size_t chunk_size, 
                               std::size_t io_block_size, 
                               std::size_t nthreads);

    // Paired binding is evaluated within one overlapped chunk.
    std::vector<AmpliconHit> search_amplicons(const std::string& path,
        const std::vector<AmpliconPrimerPair>& pairs, std::size_t threshold,
        double dg_cutoff, std::size_t amp_min, std::size_t amp_max,
        std::size_t chunk_size, std::size_t io_block_size, std::size_t nthreads);
private:
    Automaton(const std::vector<std::string>& labels,
              const std::vector<std::string>& seqs, std::size_t k,
              std::size_t threshold);
    std::vector<AmpliconHit> process_amplicons(const std::string& path,
        const std::vector<AmpliconPrimerPair>& pairs, std::size_t threshold,
        double dg_cutoff, std::size_t amp_min, std::size_t amp_max,
        std::size_t chunk_size, std::size_t io_block_size, std::size_t nthreads);
    std::vector<AmpliconHit> search_chunk_amplicons(const Task& chunk,
        const std::vector<AmpliconPrimerPair>& pairs, std::size_t owned_length,
        std::size_t threshold, double dg_cutoff, std::size_t amp_min, std::size_t amp_max);
    std::size_t search_chunk(const std::string& contig,
                             const char* data,
                             std::size_t len,
                             std::size_t owned_length,
                             std::size_t index,
                             std::size_t threshold,
                             double dg_thres,
                             std::vector<Result> &local_results);
    std::vector<Result> process_fasta_chunks(const std::string &path, 
                                             std::size_t overlap, 
                                             std::size_t threshold,
                                             double dg_thres,
                                             std::size_t chunk_size, 
                                             std::size_t io_block_size,
                                             std::size_t nthreads);
private:
    AutomatonNode *root = nullptr;
    std::vector<std::string> labels, sequences;
    std::vector<std::string> index2label, index2data;
    std::vector<std::size_t> index2pos, index2strand;
    std::unordered_map<std::string, std::size_t> label2index;
    std::size_t k;  
};
