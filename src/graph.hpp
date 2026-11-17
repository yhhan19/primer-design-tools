#pragma once
#include "utility.hpp"
#include "thal.hpp"
#include <array>

struct DimerCutoffMetrics {
    std::size_t below_cutoff = 0;
    std::size_t affected_amplicons = 0;
    std::size_t affected_primers = 0;
    double severity = 0;
};

struct DimerLandscapePoint {
    double cutoff;
    std::vector<index_t> solution;
    DimerCutoffMetrics metrics;
    std::size_t raw_affected_primers = 0;
    unsigned seed = 0;
    long time_ms = 0;
};

class KPartiteGraph {
    public:
        void set_search_seeds(const std::vector<std::vector<index_t>>& seeds);
        void enable_landscape(const std::vector<double>& cutoffs);
        void enable_adaptive_landscape(double min_cutoff, double max_cutoff);
        void finalize_landscape();
        void optimize_landscape(std::size_t points, double min_cutoff, double max_cutoff,
            bool adaptive, std::size_t restarts, std::size_t iterations, unsigned seed,
            double selection_cutoff);
        const std::vector<DimerLandscapePoint>& landscape() const { return landscape_points; }
        KPartiteGraph(const std::vector<PrimerOutput>& input, double dimer_min_c, double dimer_max_c,
                      const std::string& export_prefix = "");
        KPartiteGraph(index_t K, index_t N);
        KPartiteGraph(std::string filename);
        ~KPartiteGraph();
        void random_weights();
        index_t part(index_t i) {return i / N;}
        index_t no(index_t i) {return i % N;}
        index_t index(index_t k, index_t n) const {return k * N + n;}
        void display();
        weight_t cost(index_t *solution);
        weight_t cost(std::vector<index_t> solution);
        weight_t solve(std::size_t iter);
        std::vector<index_t> solve_fast(std::size_t iter);
        std::vector<index_t> solve_sa(std::size_t iter);
        std::vector<index_t> solve_tabu(std::size_t restarts, std::size_t iterations);
        std::vector<index_t> solve_ga(std::size_t pop_size, std::size_t generations, double mutation_rate);
        index_t get_N() { return N; }
        index_t get_K() { return K; }
        std::vector<index_t> solve_trivial();
        std::vector<index_t> solve_bottleneck(std::size_t restarts, weight_t& out_threshold);
        std::vector<index_t> solve_bottleneck_tabu(std::size_t restarts, std::size_t tenure, std::size_t max_iter, weight_t& out_threshold);
        // wrapper
        std::vector<index_t> solve_bottleneck_sa(std::size_t restarts,
                                                                std::size_t max_iter,
                                                                double      init_temp,
                                                                double      cooling,
                                                                weight_t&   out_threshold) {
            return bottleneck_search_sa(restarts, max_iter, init_temp, cooling, out_threshold);
        }
        // Maximize the (allowed_bad + 1)-th smallest selected cross-partition
        // edge. Ties at the threshold are not counted as worse. Self-dimers
        // are excluded, matching the existing bottleneck objective.
        weight_t rank_bottleneck_score(const std::vector<index_t>& solution,
                                       std::size_t allowed_bad) const;
        std::vector<index_t> solve_rank_bottleneck_sa(
            std::size_t restarts, std::size_t max_iter,
            double init_temp, double cooling, std::size_t allowed_bad,
            weight_t& out_threshold);
        std::size_t count_below_cutoff(const std::vector<index_t>& solution, weight_t cutoff) const;
        std::vector<index_t> solve_cutoff_sa(std::size_t restarts, std::size_t max_iter,
            double init_temp, double cooling, weight_t cutoff, weight_t& out_score);
        DimerCutoffMetrics dimer_cutoff_metrics(
            const std::vector<index_t>& solution, weight_t cutoff) const;
        std::vector<index_t> solve_affected_amplicons_sa(
            std::size_t restarts, std::size_t max_iter, double init_temp,
            double cooling, weight_t cutoff, weight_t& out_score);
        std::vector<index_t> solve_affected_primers_sa(
            std::size_t restarts, std::size_t max_iter, double init_temp,
            double cooling, weight_t cutoff, weight_t& out_score);
        std::vector<weight_t> selected_primer_worst_dg(
            const std::vector<index_t>& solution) const;
        weight_t primer_fraction_bound_score(
            const std::vector<index_t>& solution, std::size_t allowed_primers) const;
        std::vector<index_t> solve_primer_fraction_bound_sa(
            std::size_t restarts, std::size_t max_iter, double init_temp,
            double cooling, std::size_t allowed_primers, weight_t& out_bound);

        // Sum uses the same cross-partition edges and SA budget as rank search.
        weight_t dimer_sum_score(const std::vector<index_t>& solution) const;
        std::vector<index_t> solve_sum_sa(
            std::size_t restarts, std::size_t max_iter, double init_temp,
            double cooling, weight_t& out_score);
    private:
        std::vector<std::vector<index_t>> search_seeds;
        std::vector<DimerLandscapePoint> landscape_points;
        bool adaptive_landscape = false;
        double landscape_min = 0, landscape_max = 0;
        std::vector<double> landscape_bounds;
        std::vector<std::vector<index_t>> landscape_bound_panels;
        void observe_panel(const std::vector<index_t>& solution);
        enum class CutoffObjective { Pairs, Amplicons, Primers, PrimerBound };
        std::vector<index_t> solve_objective_sa(
            std::size_t restarts, std::size_t max_iter, double init_temp,
            double cooling, std::size_t allowed_bad, weight_t& out_score, bool sum,
            std::optional<weight_t> cutoff = std::nullopt,
            CutoffObjective cutoff_objective = CutoffObjective::Pairs);

        weight_t *vertices;
        weight_t **graph, **opt_graph, **opt_graph_2;
        index_t K, N;
        std::vector<index_t> valid_N;
        std::vector<std::vector<std::pair<index_t, index_t>>> pair_choices;
        void validate_pair_choices(const std::vector<index_t>& solution) const;
        weight_t dfs(index_t k, index_t *solution, std::size_t *count);
        void dfs2(index_t k, index_t *solution, weight_t precost, weight_t *max, std::size_t *count);
        weight_t bound(index_t k, index_t *solution);
        void preprocessing();
        weight_t random_search(std::size_t iter);
        std::vector<index_t> random_search_single_fast(std::size_t restarts);
        std::vector<index_t> beam_search_fast(std::size_t beam_width, std::size_t max_iters);
        std::vector<index_t> simulated_annealing(std::size_t restarts);
        std::vector<index_t> tabu_search(std::size_t restarts, std::size_t iterations);
        std::vector<index_t> genetic_algorithm(std::size_t pop_size, std::size_t generations, double mutation_rate);
        std::vector<index_t> bottleneck_search(std::size_t restarts, weight_t& out_threshold);
        std::vector<index_t> bottleneck_search_tabu(std::size_t restarts, std::size_t tenure, std::size_t max_iter, weight_t& out_threshold);
        weight_t             full_bottleneck(const std::vector<index_t>& sol, index_t& bn_p) const;
        weight_t             evaluate_swap(std::vector<index_t>& sol, index_t p, index_t n) const;
        // header
        std::vector<index_t> bottleneck_search_sa(std::size_t restarts, std::size_t max_iter,
                                                double init_temp, double cooling,
                                                weight_t& out_threshold, std::size_t allowed_bad = 0, bool sum = false,
                                                std::optional<weight_t> cutoff = std::nullopt,
                                                CutoffObjective cutoff_objective = CutoffObjective::Pairs);


};
