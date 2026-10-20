#include "utility.hpp"

namespace {
struct Site {
    bool known = true;
    std::size_t mismatches = 0, terminal_mismatches = 0;
    double weighted = 0, weight = 0;
};
bool acgt(char base) { return base=='A' || base=='C' || base=='G' || base=='T'; }
double mismatch_weight(std::size_t distance) {
    if(distance==0) return 20;
    if(distance==1) return 12;
    if(distance==2) return 8;
    if(distance<=4) return 4;
    if(distance<=7) return 2;
    return 1;
}
std::string site_key(const Oligo& oligo) {
    return std::to_string(oligo.start)+":"+std::to_string(oligo.length)+":"+oligo.seq;
}
bool better(const PrimerResult& a, const PrimerResult& b) {
    if(a.msa.joint_coverage != b.msa.joint_coverage)
        return a.msa.joint_coverage > b.msa.joint_coverage;
    if(a.msa.joint_3prime_match != b.msa.joint_3prime_match)
        return a.msa.joint_3prime_match > b.msa.joint_3prime_match;
    return a.msa.weighted_mismatch < b.msa.weighted_mismatch;
}
}

void score_candidate_pairs_msa(PrimerOutput& output,
    const std::vector<std::string>& sequences, const Args& args) {
    if(sequences.empty()) throw std::invalid_argument("Pair MSA scoring requires aligned sequences");
    if(args.pair_3prime_window==0) throw std::invalid_argument("Pair MSA 3-prime window must be positive");
    // Cache each site's per-sequence assessment once, then join partner
    // assessments by MSA sequence index. Do not multiply marginal coverage.
    std::unordered_map<std::string, std::vector<Site>> left_cache, right_cache;
    auto assess = [&](const Oligo& oligo, bool right) -> const std::vector<Site>& {
        auto& cache=right ? right_cache : left_cache;
        const auto key=site_key(oligo);
        auto found=cache.find(key);
        if(found!=cache.end()) return found->second;
        std::vector<Site> sites(sequences.size());
        for(std::size_t i=0;i<sequences.size();++i) {
            auto& site=sites[i];
            if(oligo.start<0 || oligo.length<=0 || oligo.seq.size()!=std::size_t(oligo.length) ||
               std::size_t(oligo.start)+oligo.length>sequences[i].size()) {
                site.known=false; continue;
            }
            for(std::size_t j=0;j<std::size_t(oligo.length);++j) {
                const char a=oligo.seq[j], b=sequences[i][oligo.start+j];
                if(!acgt(a) || !acgt(b)) {site.known=false; continue;}
                const auto distance=right ? j : oligo.length-1-j;
                const auto weight=mismatch_weight(distance);
                site.weight += weight;
                if(a!=b) {
                    ++site.mismatches; site.weighted += weight;
                    if(distance<args.pair_3prime_window) ++site.terminal_mismatches;
                }
            }
        }
        return cache.emplace(key,std::move(sites)).first->second;
    };
    for(auto& pair:output.pairs) {
        const auto& left=assess(pair.left,false);
        const auto& right=assess(pair.right,true);
        std::size_t passing=0, perfect=0, terminal=0, unknown=0;
        double weighted=0, weight=0;
        for(std::size_t i=0;i<sequences.size();++i) {
            const auto& l=left[i]; const auto& r=right[i];
            if(!l.known || !r.known) {++unknown;continue;}
            if(l.mismatches<=args.pair_max_mismatches && r.mismatches<=args.pair_max_mismatches &&
               l.terminal_mismatches<=args.pair_max_3prime_mismatches &&
               r.terminal_mismatches<=args.pair_max_3prime_mismatches) ++passing;
            if(l.mismatches==0 && r.mismatches==0) ++perfect;
            if(l.terminal_mismatches==0 && r.terminal_mismatches==0) ++terminal;
            weighted += l.weighted+r.weighted;
            weight += l.weight+r.weight;
        }
        const double n=sequences.size();
        pair.msa.joint_coverage=passing/n;
        pair.msa.joint_perfect_match=perfect/n;
        pair.msa.joint_3prime_match=terminal/n;
        pair.msa.unknown_fraction=unknown/n;
        pair.msa.weighted_mismatch=weight>0 ? weighted/weight : 1;
        pair.msa.sequence_count=sequences.size();
    }
    std::stable_sort(output.pairs.begin(),output.pairs.end(),better);
}

std::size_t retain_candidate_pairs_msa(PrimerOutput& output, const Args& args) {
    if(!std::isfinite(args.pair_keep_percent) || args.pair_keep_percent<0 ||
       args.pair_keep_percent>100 || args.pair_min_return==0)
        throw std::invalid_argument("Invalid pair MSA retention settings");
    if(output.pairs.empty()) return 0;
    for(const auto& pair:output.pairs)
        if(pair.msa.sequence_count==0) throw std::invalid_argument("Pair MSA scores must be computed before retention");
    const auto percent_count=static_cast<std::size_t>(std::ceil(
        static_cast<long double>(output.pairs.size()) * args.pair_keep_percent / 100.0L));
    std::size_t keep=std::min(output.pairs.size(),std::max(args.pair_min_return,percent_count));
    // Scores are sorted by score_candidate_pairs_msa. Keep complete boundary ties.
    while(keep<output.pairs.size() && !better(output.pairs[keep-1],output.pairs[keep])) ++keep;
    output.pairs.resize(keep);
    auto remove_orphans=[&](std::vector<Oligo>& oligos, bool right) {
        std::set<std::string> used;
        for(const auto& pair:output.pairs) used.insert(site_key(right ? pair.right : pair.left));
        oligos.erase(std::remove_if(oligos.begin(),oligos.end(),[&](const Oligo& o) {
            return !used.count(site_key(o));
        }),oligos.end());
    };
    remove_orphans(output.left_oligos,false);
    remove_orphans(output.right_oligos,true);
    return keep;
}
