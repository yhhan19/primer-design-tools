#include "utility.hpp"

bool same_oligo(const Oligo& a, const Oligo& b) {
    return a.start == b.start && a.length == b.length && a.seq == b.seq;
}

void retain_surviving_pairs(PrimerOutput& output) {
    auto present = [](const Oligo& oligo, const std::vector<Oligo>& candidates) {
        return std::any_of(candidates.begin(), candidates.end(),
            [&](const Oligo& candidate) { return same_oligo(oligo, candidate); });
    };
    output.pairs.erase(std::remove_if(output.pairs.begin(), output.pairs.end(),
        [&](const PrimerResult& pair) {
            return !present(pair.left, output.left_oligos) ||
                   !present(pair.right, output.right_oligos);
        }), output.pairs.end());
}

PairingStats build_candidate_pairs(PrimerOutput& output, const std::string& target, const Args& args) {
    PairingStats stats;
    stats.combinations = output.left_oligos.size() * output.right_oligos.size();
    output.pairs.clear();
    output.pairs_required = true;
    using Settings = std::unique_ptr<p3_global_settings, decltype(&p3_destroy_global_settings)>;
    using Sequence = std::unique_ptr<seq_args, decltype(&destroy_seq_args)>;
    using ResultPtr = std::unique_ptr<p3retval, decltype(&destroy_p3retval)>;
    Settings settings(p3_create_global_settings(), p3_destroy_global_settings);
    if (!settings) throw std::runtime_error("Cannot allocate Primer3 pairing settings");
    auto* pa = settings.get();
    p3_set_gs_primer_opt_size(pa, args.p3_opt_size);
    p3_set_gs_primer_min_size(pa, args.p3_min_size);
    p3_set_gs_primer_max_size(pa, args.p3_max_size);
    p3_set_gs_primer_opt_tm(pa, args.p3_opt_tm);
    p3_set_gs_primer_min_tm(pa, args.p3_min_tm);
    p3_set_gs_primer_max_tm(pa, args.p3_max_tm);
    p3_set_gs_primer_min_gc(pa, args.p3_min_gc);
    p3_set_gs_primer_max_gc(pa, args.p3_max_gc);
    pa->p_args.salt_conc = args.mv;
    pa->p_args.divalent_conc = args.dv;
    pa->p_args.dntp_conc = args.dntp;
    pa->p_args.dna_conc = args.dna_conc;
    pa->primer_task = generic;
    pa->pick_left_primer = 1;
    pa->pick_right_primer = 1;
    pa->pick_internal_oligo = 0;
    pa->pick_anyway = 0;
    pa->num_return = 1;
    pa->num_intervals = 1;
    // PDR layout limits do not filter candidate combinations. Let Primer3
    // assess pairs over its legal product interval for this target.
    if (target.size() > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Pairing target exceeds Primer3 coordinate range");
    pa->pr_min[0] = args.p3_max_size;
    pa->pr_max[0] = static_cast<int>(target.size());
    for (const auto& left : output.left_oligos) {
        for (const auto& right : output.right_oligos) {
            if (left.start < 0 || right.start < 0 || left.length <= 0 || right.length <= 0)
                continue;
            const auto left_end = std::size_t(left.start) + left.length;
            const auto right_end = std::size_t(right.start) + right.length;
            if (left_end > target.size() || right_end > target.size())
                throw std::invalid_argument("Candidate binding interval is outside the target");
            if (right_end <= std::size_t(left.start))
                throw std::invalid_argument("Candidate right endpoint precedes left start");
            const auto product_size = right_end - std::size_t(left.start);
            // Oligo.seq stores the reference-oriented segment for both sides.
            // Primer3's explicit right primer input must instead be 5'->3'.
            if (target.compare(left.start, left.length, left.seq) != 0 ||
                target.compare(right.start, right.length, right.seq) != 0)
                throw std::runtime_error("Candidate sequence does not match its target interval");
            ++stats.checked;
            Sequence sequence(create_seq_arg(), destroy_seq_args);
            if (!sequence) throw std::runtime_error("Cannot allocate Primer3 pairing sequence");
            auto* sa = sequence.get();
            sa->sequence = strdup(target.c_str());
            sa->sequence_name = strdup("candidate_pair");
            p3_set_sa_left_input(sa, left.seq.c_str());
            const std::string right_primer = revcomp(right.seq);
            p3_set_sa_right_input(sa, right_primer.c_str());
            sa->force_left_start = left.start;
            sa->force_left_end = left.start + left.length - 1;
            sa->force_right_start = right.start + right.length - 1;
            sa->force_right_end = right.start;
            ResultPtr result(choose_primers(pa, sa), destroy_p3retval);
            if (!result) throw std::runtime_error("Primer3 pairing returned NULL");
            if (result->glob_err.data && result->glob_err.data[0])
                throw std::runtime_error(std::string("Primer3 pairing: ") + result->glob_err.data);
            if (result->per_sequence_err.data && result->per_sequence_err.data[0])
                throw std::runtime_error(std::string("Primer3 pairing: ") + result->per_sequence_err.data);
            if (result->best_pairs.num_pairs == 0) continue;
            const auto& validated = result->best_pairs.pairs[0];
            if (validated.left->start != left.start ||
                validated.left->length != left.length ||
                validated.right->start != right.start + right.length - 1 ||
                validated.right->length != right.length ||
                validated.product_size != static_cast<int>(product_size))
                throw std::runtime_error("Primer3 returned a different candidate pair");
            PrimerResult pair;
            pair.left = left;
            pair.right = right;
            pair.product_size = product_size;
            pair.pdr_left = output.pdr_left;
            pair.pdr_right = output.pdr_right;
            pair.segment_id = output.segment_id;
            output.pairs.push_back(std::move(pair));
        }
    }
    stats.accepted = output.pairs.size();
    return stats;
}
