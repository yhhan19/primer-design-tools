#include "automaton.hpp"
#include <tuple>
#include <array>

std::vector<AmpliconHit> Automaton::search_chunk_amplicons(const Task& chunk,
    const std::vector<AmpliconPrimerPair>& pairs, std::size_t owned_length,
    std::size_t threshold, double dg_cutoff, std::size_t amp_min, std::size_t amp_max) {
    std::vector<Result> hits;
    search_chunk(chunk.contig, chunk.seq.data(), chunk.seq.size(), 0, threshold, dg_cutoff, hits);
    // Single-hit search returns positions local to this chunk. Pair locally,
    // then translate both binding sites and the product to global coordinates.
    std::unordered_map<std::string, std::array<std::vector<const Result*>,2>> by_label;
    for(const auto& hit:hits) by_label[hit.label][hit.flag].push_back(&hit);
    for(auto& entry:by_label) for(auto& strand:entry.second)
        std::sort(strand.begin(),strand.end(),[](const Result* a,const Result* b){return a->ref_index<b->ref_index;});
    std::vector<AmpliconHit> products;
    for(std::size_t pair_id=0;pair_id<pairs.size();++pair_id) {
        const auto& pair=pairs[pair_id];
        auto first=by_label.find(pair.first_label), second=by_label.find(pair.second_label);
        if(first==by_label.end() || second==by_label.end()) continue;
        auto join=[&](const std::vector<const Result*>& forward,
                      const std::vector<const Result*>& reverse) {
            for(const auto* left:forward) {
                if(left->ref_index>=owned_length) continue; // Exactly one owning chunk.
                const auto left_end=left->ref_index+left->data.size();
                auto begin=std::lower_bound(reverse.begin(),reverse.end(),left_end,
                    [](const Result* hit,std::size_t pos){return hit->ref_index<pos;});
                for(auto it=begin;it!=reverse.end();++it) {
                    const auto* right=*it;
                    if(right->ref_index-left->ref_index>=amp_max) break;
                    const auto end=right->ref_index+right->data.size();
                    const auto size=end-left->ref_index;
                    if(size<amp_min || size>amp_max) continue;
                    AmpliconHit product;
                    product.pair_index=pair_id;product.contig=chunk.contig;
                    product.start=chunk.start+left->ref_index;product.end=chunk.start+end;
                    product.forward_hit=*left;product.reverse_hit=*right;
                    product.forward_hit.ref_index+=chunk.start;
                    product.reverse_hit.ref_index+=chunk.start;
                    products.push_back(std::move(product));
                }
            }
        };
        // Stored right candidates use reference-oriented sequence: their
        // physical extension direction is the opposite of a left candidate's.
        join(first->second[pair.first_is_right ? 1 : 0], second->second[pair.second_is_right ? 0 : 1]);
        // Also search the reversed genomic arrangement of the same pair.
        join(second->second[pair.second_is_right ? 1 : 0], first->second[pair.first_is_right ? 0 : 1]);
    }
    return products;
}

std::vector<AmpliconHit> Automaton::search_amplicons(const std::string& path,
    const std::vector<AmpliconPrimerPair>& pairs, std::size_t threshold,
    double dg_cutoff, std::size_t amp_min, std::size_t amp_max,
    std::size_t chunk_size, std::size_t io_block_size, std::size_t nthreads) {
    if(!amp_min || amp_min>amp_max || !chunk_size || !io_block_size || !nthreads ||
       nthreads>std::numeric_limits<std::size_t>::max()/2 || !std::isfinite(dg_cutoff))
        throw std::invalid_argument("Invalid paired off-target search settings");
    const auto overlap=amp_max-1;
    if(overlap>std::numeric_limits<std::size_t>::max()-chunk_size)
        throw std::invalid_argument("Paired off-target chunk size overflow");
    for(const auto& pair:pairs)
        if(!label2index.count(pair.first_label) || !label2index.count(pair.second_label))
            throw std::invalid_argument("Paired off-target search refers to an unknown primer label");
    std::ifstream input(path,std::ios::binary);
    if(!input) throw std::runtime_error("Cannot open off-target reference: "+path);
    std::vector<char> io(io_block_size);
    TaskQueue queue(std::max<std::size_t>(2,nthreads*2));
    std::vector<std::vector<AmpliconHit>> buckets(nthreads);
    std::vector<std::thread> workers;
    std::mutex error_mutex;
    std::exception_ptr error;
    std::atomic<bool> failed{false};
    std::atomic<std::size_t> processed{0};
    auto capture_error=[&]() {
        std::lock_guard<std::mutex> lock(error_mutex);
        if(!error) error=std::current_exception();
        failed=true;
    };
    try {
        for(std::size_t t=0;t<nthreads;++t) workers.emplace_back([&,t]() {
            Task job;
            while(queue.pop(job)) {
                if(failed) continue;
                try {
                    auto products=search_chunk_amplicons(job,pairs,
                        std::min(chunk_size,job.seq.size()),threshold,dg_cutoff,amp_min,amp_max);
                    auto& bucket=buckets[t];
                    bucket.insert(bucket.end(),std::make_move_iterator(products.begin()),std::make_move_iterator(products.end()));
                    ++processed;
                } catch(...) {capture_error();}
            }
        });
    } catch(...) {
        queue.close();for(auto& worker:workers) worker.join();throw;
    }
    std::size_t contigs=0,chunks=0,offset=0;
    std::string contig,header,buffer;
    std::set<std::string> contig_names;
    bool in_header=false,have_contig=false,saw_sequence=false;
    auto last_progress=std::chrono::steady_clock::now();
    auto flush=[&]() {
        if(have_contig && !buffer.empty()) {
            queue.push(Task{contig,std::move(buffer),offset,true});++chunks;
            buffer.clear();
        }
    };
    try {
        std::cout << "Paired off-target chunk search: " << pairs.size() << " retained pairs, products "
                  << amp_min << "-" << amp_max << " bp, core " << chunk_size
                  << " bp + " << overlap << " bp overlap, " << nthreads << " threads\n";
        while(input && !failed) {
            input.read(io.data(),io.size());
            const auto count=input.gcount();
            for(std::streamsize i=0;i<count && !failed;++i) {
                const unsigned char raw=static_cast<unsigned char>(io[i]);
                if(in_header) {
                    if(raw=='\n' || raw=='\r') {
                        contig=header.substr(0,header.find_first_of(" \t"));
                        if(contig.empty() || !contig_names.insert(contig).second)
                            throw std::runtime_error("Off-target FASTA requires unique nonempty contig names");
                        header.clear();in_header=false;have_contig=true;offset=0;++contigs;
                    } else header.push_back(raw);
                    continue;
                }
                if(raw=='>') {flush();header.clear();in_header=true;have_contig=false;continue;}
                if(std::isspace(raw)) continue;
                if(!have_contig) throw std::runtime_error("Off-target reference is not valid FASTA");
                const char base=std::toupper(raw);
                buffer.push_back(base=='A'||base=='C'||base=='G'||base=='T' ? base : 'N');
                saw_sequence=true;
                if(buffer.size()>=chunk_size+overlap) {
                    queue.push(Task{contig,buffer,offset,false});++chunks;
                    buffer.erase(0,chunk_size);offset+=chunk_size;
                }
                if ((i & 65535) == 0) {
                    const auto now=std::chrono::steady_clock::now();
                    if(now-last_progress>=std::chrono::seconds(10)) {
                        std::cout << "Paired reference search: " << contigs << " contigs, "
                                  << processed.load() << "/" << chunks << " chunks processed\n";
                        last_progress=now;
                    }
                }
            }
        }
        if(input.bad()) throw std::runtime_error("Error reading off-target reference");
        if(!saw_sequence) throw std::runtime_error("Off-target reference contains no sequence");
        flush();
    } catch(...) {capture_error();}
    queue.close();for(auto& worker:workers) worker.join();
    if(error) std::rethrow_exception(error);
    std::vector<AmpliconHit> results;
    for(auto& bucket:buckets) results.insert(results.end(),std::make_move_iterator(bucket.begin()),std::make_move_iterator(bucket.end()));
    auto key=[](const AmpliconHit& h) {
        return std::make_tuple(h.contig,h.start,h.end,h.pair_index,
            h.forward_hit.label,h.forward_hit.flag,h.reverse_hit.label,h.reverse_hit.flag);
    };
    std::sort(results.begin(),results.end(),[&](const AmpliconHit& a,const AmpliconHit& b){return key(a)<key(b);});
    results.erase(std::unique(results.begin(),results.end(),[&](const AmpliconHit& a,const AmpliconHit& b){return key(a)==key(b);}),results.end());
    std::cout << "Paired reference search completed: " << contigs << " contigs, " << chunks
              << " chunks, " << results.size() << " potential products\n";
    return results;
}
