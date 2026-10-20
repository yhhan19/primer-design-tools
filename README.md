# Δ-PRO: Primer Design and Off-Target Analysis Tool

Δ-PRO is a comprehensive C++ tool for multiplex PCR primer design optimization and tiling that implements a four-stage pipeline: PDR (Primer Design Region) optimization to identify high-quality, conserved regions with sufficient primer candidates; primer selection using the established Primer3 algorithm to generate high-quality primer pairs with proper thermodynamic constraints; off-target screening against reference databases using k-mer based searching to filter primers with potential off-target binding; and dimer minimization using multiple optimization algorithms to select primer sets with minimal primer-dimer formation potential, ensuring both specificity and compatibility for successful multiplex PCR amplification and tiling coverage.

## Quick Start

The fastest way to run Δ-PRO is to clone the repository, build the `dpro` executable with a local Primer3 installation, and run the pipeline on a FASTA input file.

```bash
# Clone the repository
git clone https://github.com/yhhan19/primer-design-tools.git
cd primer-design-tools

# Build Δ-PRO
make PRIMER3=/path/to/primer3

# Run the pipeline using example data
./bin/dpro \
  -i ./data/sequence.fasta \
  -o ./data/results
```

**Note:** Replace `/path/to/primer3` with the actual path to your [Primer3](https://github.com/primer3-org/primer3/) installation directory. For example, if Primer3 is located in the repository directory as `./primer3`, build Δ-PRO with:

```bash
make PRIMER3=./primer3
```

The final primer solution will be written to:

```bash
./data/results.bed
```

## Installation

### Prerequisites
- **C++17 compatible compiler** (GCC 7+, Clang 5+)
- **Primer3 library** (https://github.com/primer3-org/primer3)

### Build
```bash
git clone https://github.com/yhhan19/primer-design-tools.git
cd primer-design-tools
make PRIMER3=/path/to/primer3
```
**Note:** Replace `/path/to/primer3` with the actual path to the Primer3 installation directory.
The executable will be generated as `bin/dpro`.

## Usage

```bash
./bin/dpro -i <input> -o <output> [options]
```

### Core Arguments
```bash
-i <file>    Input sequence file (required)
-o <file>    Output directory (required)
-r <file>    Reference genome file (optional)
-h, --help   Show help message
```

### Other Parameters

#### Search Parameters
```bash
-k <int>       Kmer length for indexing (default: 15)
-H <int>       Hamming distance threshold (default: 3)  
-H2 <int>      Alternative Hamming threshold (default: 3)
-t <int>       Number of threads (default: 8)
-C <int>       Processing chunk size (default: 100000)
-B <int>       I/O block size (default: 8192)
```

#### Optimization Parameters
```bash
-Ln <int>     Maximum amplicon length (default: 420)
-Lx <int>     Minimum amplicon length (default: 252)  
-Lp <int>     Primer design region length (default: 40)
--pdr-search <fixed|midpoint>
              PDR layout search (default: fixed)
-Ux <double>  Maximum risk threshold (default: 10000.0)
-Un <double>  Minimum risk threshold (default: 0.1)
-I <int>      Optimization iterations (default: 1000)
-S <int>      Random seed (default: 42)
```

`--pdr-search fixed` preserves the existing efficient search and gives every
amplicon the `-Ln` length. `--pdr-search midpoint` permits each amplicon to
vary from `-Lx` through `-Ln` and applies the midpoint/non-overlap recurrence
between adjacent amplicons. Its dynamic program uses a two-dimensional range
minimum query indexed by predecessor right coordinate and midpoint sum, giving
`O(log^2 n)` predecessor updates and queries instead of scanning the full
two-dimensional predecessor region.

#### Thermodynamic Parameters  
```bash
-G <double>     dG threshold for filtering (default: -20000.0)
--mv <double>   Monovalent concentration mM (default: 50.0)
--dv <double>   Divalent concentration mM (default: 1.5) 
--dntp <double> dNTP concentration mM (default: 0.6)
--dna <double>  DNA concentration nM (default: 50.0)
--temp <double> Temperature °C (default: 60.0)
```

## Examples

### Basic Pipeline Analysis
```bash
# Complete pipeline with reference genome screening
./bin/dpro \
  -i sequences.fasta \
  -o results/ \
  -r reference_genome.fna

# Pipeline without off-target screening
./bin/dpro \
  -i sequences.fasta \
  -o results/

# Custom parameters for stringent filtering
./bin/dpro \
  -i sequences.fasta \
  -o results/ \
  -r genome.fna \
  -G -25000.0 \
  -t 16 \
  -H 2
```

### Advanced Configuration
```bash
# High-throughput analysis
./bin/dpro \
  -i large_dataset.fasta \
  -o analysis_results/ \
  -r combined_genomes.fna \
  -t 32 \

# Custom thermodynamic conditions
./bin/dpro \
  -i primers.fasta \
  -o filtered_results/ \
  -r host_genome.fna \
  --temp 65.0 \
  --mv 100.0 \
  --dna 25.0
```

## Input Format
- **Input**: FASTA format with target sequences
- **Reference**: FASTA format reference genomes (optional)

## Output Structure
### Pipeline Stages
1. **PDR (Primer Design Region) Optimization** - Identifies high-quality, conserved regions with sufficient primer candidates and achieves optimal tiling coverage through constrained optimization
2. **Primer Selection** - Uses Primer3 algorithm to generate high-quality primer pairs with proper thermodynamic constraints within optimized regions
3. **Off-target Screening** - Reference genome analysis using k-mer based searching with thermodynamic filtering to remove primers with potential off-target binding
4. **Dimer Minimization** - Employs multiple optimization algorithms (Random Search, Simulated Annealing, Tabu Search, Genetic Algorithm) to select primer sets with minimal primer-dimer formation potential

### Output Structure

The pipeline generates comprehensive analysis reports:

#### PDR Optimization Results
```
PDRs: 20
uncovered: 3 (front), 308 (rear)
loss: 109.869

Search for optimal PDR...
--------------------------------------------
    iter          u'       loss'        loss
--------------------------------------------
       1        5000        9000     257.973
       2        2500        4500     257.973
       ......   ......      ......    ......
      17     52.5665     109.869     109.869
--------------------------------------------
Search completed!
```

#### Primer Selection Summary
```
Starting primer selection for 10 PDRs
Product size: [252, 460] bp  |  Template length: 3063 bp
------------------------------------------------------------
     Index                PDRs      Left     Right     Pairs
------------------------------------------------------------
         0            [3, 383]        46        73         5
         1          [252, 632]        80        19         5
         2          [585, 965]        47        87         5
         3         [800, 1180]        50        64         5
         4        [1122, 1502]        49        20         5
    ......              ......    ......    ......    ......
------------------------------------------------------------
```

#### Off-target Screening Results
```
Searching ./scripts/ref/livestock_combined.fna (threads: 8)
----------------------------------------------------------
   Contigs      Chunks    Candidates   Time(s)  Rate(ch/s)
----------------------------------------------------------
         6          13           269        11       14.1M
        29          41          2929        35       44.7M
    ......       ......       ......    ......      ......
      2358        2387         10943       116       47.0M
----------------------------------------------------------
Found 16756 candidates in reference
Total primers to remove: 5057/5057 (dG <= -15000)

Summary
--------------------------------------------------------------------------------
PDR pair  Original (P/L/R)   Removed (P/L/R)     Final (P/L/R)         Threshold
--------------------------------------------------------------------------------
       0           5/46/73           5/16/27           0/30/46    -15000.000000
       1           5/80/19           5/27/14            0/53/5    -15000.000000
  ......            ......            ......            ......           ......
--------------------------------------------------------------------------------
   TOTAL               930               411               519
--------------------------------------------------------------------------------
```

#### Dimer Minimization Results
```
--------------------------------------------------
Algorithm                 Cost  Time(ms)    Winner
--------------------------------------------------
Random Search       -3223.0680      2485       ✓
SA                  -3223.0680     65825
Tabu Search         -3223.0680      5422
Genetic             -3223.0680     32909
--------------------------------------------------
 Winner: Random Search  cost=-3223.0680  (2485 ms)
```

#### Final Solution Primers
```
Solution Primers
-------------------------------------------------------------------------------------------------------------------
   idx             PDR                    left seq     Tm    GC%                   right seq     Tm    GC%   size
-------------------------------------------------------------------------------------------------------------------
     0        [3, 383]        TTAAAATAGTGTCGCCGACG   58.2   45.0      GGATGTTCAGCCTCTAAAGGTT   59.9   45.5    419
     1      [252, 632]     CGAGTGTTAGATATCACACTGAG   58.0   43.5   TATAATAGTGAACGCTGAAAGGAGG   60.1   40.0    411
     2      [585, 965]      GTATTGAAATGTTTCACGCAGC   58.9   40.9      GATACGGAGGAGGAACCTAAGG   60.9   54.5    403
     ......     ......                      ...... ...... ......                      ...... ...... ...... ......
-------------------------------------------------------------------------------------------------------------------
```

## License

This project is licensed under the GNU General Public License v3.0 or later.
See the [LICENSE](LICENSE) file for details.

## Citation
If you use Δ-PRO in your research, please cite: 
Scaling Variant-Aware Multiplex Primer Design 
Yunheng Han, Christina Boucher 
bioRxiv 2026.02.03.703607
doi: https://doi.org/10.64898/2026.02.03.703607
```
@article{han2026scaling,
  title   = {Scaling Variant-Aware Multiplex Primer Design},
  author  = {Han, Yunheng and Boucher, Christina},
  journal = {bioRxiv},
  year    = {2026},
  doi     = {10.64898/2026.02.03.703607},
  url     = {https://doi.org/10.64898/2026.02.03.703607}
}
```

## Dimer annealing-window scoring

Add `--dimer-temp-min 58 --dimer-temp-max 62` to score graph edges by the
most negative predicted dimer free energy in that inclusive Celsius window.
These are example settings: choose the window for the intended PCR protocol.
Both options are required together; equal endpoints select one temperature.
Omitting both preserves the existing `--temp` behavior (37 C by default).
Saved argument files use `dimer_temp_min` and `dimer_temp_max`.

For the structure returned by libthal, the weight is
`min(dH - (Tmin + 273.15)*dS, dH - (Tmax + 273.15)*dS)` in cal/mol.
This endpoint calculation is exact for constant dH/dS; for negative dS the
minimum is at Tmin. More negative weights remain worse, so the existing
bottleneck maximization retains its direction. Only dimer graph edge scoring
uses this range; off-target/candidate calculations still use `--temp`.

This is a fixed-structure thermodynamic score. The bundled libthal selects
dimer structures using a reference temperature; this change does not perform
a temperature-dependent structure search, predict dimer yield, add self-dimer
edges, or establish that intended binding is adequate throughout the window.

## Unweighted K-partite graph export

Each dimer pool automatically exports its candidate graph before calculating
that pool's dimer weights. With `-o results`, files are:

- `results.pool1.kpartite.nodes.tsv` (and pool2): real candidate nodes,
  partition IDs, pool-local amplicon/candidate numbers, segment ID, L/R side,
  sequence, and original stored coordinates.
- `results.pool1.kpartite.edges.tsv`: each undirected cross-partition edge
  once, with `NA` for the not-yet-computed weight (not zero).
- `results.pool1.kpartite.meta.tsv`: counts, ID conventions, temperatures,
  and unweighted status. Empty partitions remain represented by the partition
  count and ID convention: partition 2*a is left, 2*a+1 right.

Node IDs match internal graph IDs (`partition_id * node_id_stride + candidate_index`),
so gaps are intentional when partitions have unequal candidate counts.
No padded nodes, self-edges, or within-partition edges are exported.
The same output prefix overwrites these reports on subsequent runs. Explicit
edge files grow quadratically with candidate count. Files are closed and checked
before scoring starts; export failures stop the stage. These TSVs are an
unweighted interchange format, not the legacy weighted graph-loader format.

After a pool finishes dimer scoring, its `.edges.tsv` is rewritten with the
calculated range-minimum dG weights in cal/mol (full double precision), and
`.meta.tsv` status becomes `weighted`. Before scoring completes, weights are
`NA` and status is `unweighted`. The node IDs and edge order do not change.
Check metadata status before treating an export as a completed weighted graph.

### Fixed-rank dimer optimization

Dimer selection uses simulated annealing to maximize the `(k+1)`-th smallest
selected interaction dG in each pool, with `k = floor(M * percent / 100)`
and `M = n * (n-1) / 2` for `n` selected primers. Set the percentage with
`--dimer-bad-percent 1` (default: 1%; finite values in `[0,100)`).
Use `--dimer-bad-percent 0` for the original worst-edge objective.
The saved `.args` key is `dimer_bad_percent`. Scores use cal/mol; more
negative dG is worse. Only values strictly below the threshold count as
worse, so ties can leave fewer than k interactions below it. Self-dimers
are excluded from this objective. Each pool reports its edge count, k,
and optimized threshold. This is a heuristic search, not a guarantee of
a globally optimal panel.

### Comparing dimer objectives

Use `--dimer-methods sum,min-max,top-k-bound` to run multiple objectives.
The default remains `top-k-bound`. Each method uses simulated annealing
with the same candidate graph, initial temperature, cooling schedule,
4000 restarts, `5 * -I` proposals per restart, and seed reset per method.
`sum` maximizes summed selected cross-primer dG; `min-max` maximizes the
worst selected dG; `top-k-bound` maximizes the (k+1)-th smallest dG using
`--dimer-bad-percent`. Self-dimers are excluded from all three objectives.

The first listed method supplies the final panel; override this with
`--dimer-select min-max` (must be among the listed methods). Scores from
different objectives are not treated as directly comparable winners.
Every result is evaluated under all three metrics, where higher is better,
and timing excludes graph construction and report writing. Identical search
settings do not guarantee identical convergence across objectives.

Per pool, `<output>.poolN.dimer-methods.tsv` records the comparison, and
`<output>.poolN.<method>.primers.tsv` preserves each alternative's selected
primer sequences and indices. Target indices in these files are pool-local.
The saved configuration keys are `dimer_methods` and `dimer_select`.
The prior same-selected-primers dimer score comparison remains disabled.

### Current run defaults

`./bin/dpro -i ./data/sequence.fasta -o ./data/results` now uses
`--p3-weighted-keep-percent 20 --pdr-search midpoint --p3-num-return 5`
and `--dimer-temp-min 58 --dimer-temp-max 62` plus
`--dimer-methods sum,min-max,top-k-bound --dimer-bad-percent 1`
with `top-k-bound` supplying the final panel. All options remain overridable.
Without an explicit `--dimer-select`, top-k-bound is preferred when listed;
otherwise the first listed method is selected. Explicit selections must be
among the requested methods. Dimer temperatures default to 58-62 C even
when `--temp` is changed; use equal dimer min/max values for a single
dimer scoring temperature. Each endpoint can be overridden independently.
These defaults supersede earlier default descriptions in this document.

### Parameter reporting

`./bin/dpro --help` documents every accepted CLI option and reads defaults
from the same Args settings as execution. The startup log reports all
effective settings, configuration files, optional reference/risk inputs,
selected pipeline mode, mismatch mode, PDR forward/reverse limits, and
fixed/derived dimer SA settings. Use `--args <output>.args` to replay saved
settings; later CLI options or configurations override earlier values.
The optional `-a` risk file switches execution to PDR only. `-H` uses the
standard mismatch count; `-H2` uses the existing adjacent-mismatch count.
Saved configurations preserve this with `threshold_mode = standard` or
`adjacent`, alongside the user-level `threshold` value. The off-target
`-G` default is -15000 cal/mol. Concentrations are mM except DNA in nM.

### Candidate pairing before downstream stages

After per-side variant ranking and retention, every retained left/right
combination is checked for inward-facing non-overlapping target geometry,
actual product length within -Lx..-Ln, and Primer3 pair acceptance using
the configured individual-primer settings and Primer3 default pair criteria.
Exact candidate positions are forced; pair checking cannot substitute a
different binding site. Reference-oriented right segments are reverse
complemented when supplied as explicit Primer3 right-primer inputs.

Validated pairs preserve candidate variant metrics and target metadata.
They are exported to `<output>.segmentN.pdrN.pairs.tsv` (zero-based indices).
The right sequence column retains the existing reference-oriented convention.
Single-primer off-target screening remains in use; it removes affected
validated combinations and never rebuilds the unrestricted Cartesian product.
If no valid pair remains, the pipeline reports infeasibility explicitly.
All three CLI dimer objectives initialize and mutate whole validated pairs,
while scoring the same primer-level interaction edges as before. The graph
export includes `.allowed_pairs.tsv` describing the additional constraints.
This pairing change does not implement the proposed paired off-target scan.

### Variation filtering temporarily disabled

All individual Primer3 candidates now reach explicit pair validation.
Variation metrics remain available for ranking and reporting, but neither
the adaptive variation-score break nor the weighted keep percentage removes
candidates. `--p3-weighted-keep-percent` remains accepted and saved for
configuration compatibility; help and logs mark it inactive. Primer3 pair
checks, product-length constraints and off-target filtering still apply.
This may increase pairing and dimer-search runtime with larger candidate sets.

### Geometry prefilter removed

Candidate pairing no longer rejects combinations for overlap or product
length outside -Lx..-Ln. Those limits control PDR layout only. Every valid
candidate combination is submitted to Primer3 using the broad legal product
interval from PRIMER_MAX_SIZE to the target length. Primer3 may still reject
combinations under its own pair criteria. Basic coordinate checks remain to
prevent invalid target access. The pair summary now reports P3 checked in
place of Geom rejected. This supersedes earlier pairing-length descriptions.

### MSA filtering of validated pairs

After Primer3 checks, pairs are scored against aligned MSA rows and ranked
by joint coverage (descending), joint exact 3-prime match (descending),
then combined weighted mismatch (ascending). Joint coverage requires both
primers to pass on the same row; it is not a product of marginal coverage.
Defaults: `--pair-max-mismatches 2 --pair-max-3prime-mismatches 0` with
`--pair-3prime-window 5`. Right-primer 3-prime coordinates account for the
existing reference-oriented stored sequence. All rates use the total MSA
row count; a gap, ambiguity or unavailable base in either full binding site
makes the pair unknown for that row, which contributes to neither coverage
nor match success. Weighted mismatch uses only rows with both known sites;
all-unknown pairs receive weighted score 1.

`--pair-keep-percent 20 --pair-min-return 5` retains the larger of the rounded
up percentage and minimum, capped by available valid pairs, including all
boundary ties. A minimum safeguards candidate count, not a required biological
coverage threshold. Use 100 percent to retain every valid pair. Individual
variation-based filtering remains disabled. Primers without retained pair
partners are removed from downstream candidate lists. Each target writes
`.pair-msa.tsv` with all ranked scores and keep/drop flags, while `.pairs.tsv`
contains retained pairs. The summary separates P3 valid from MSA kept.

### Paired off-target chunk search

Off-target filtering now defaults to `--off-target-mode paired`, with
`--off-target-amp-min 1 --off-target-amp-max 1000` in bp. It screens the
retained candidate pairs against `-r`, reusing the existing seed/mismatch
and thermodynamic hit checks. A pair is rejected only when BOTH hits have
dG strictly below `-G`, on the same contig, with inward-facing non-overlapping
binding sites and a full product length in the configured inclusive range.
Both genomic arrangements of each candidate pair are searched.

The new `search_chunk_amplicons` pairs hits within ONE chunk; cores have
max-product-minus-one bases of right overlap. Only the chunk owning the
leftmost product coordinate reports it. Coordinates are global zero-based
half-open intervals. Results are deterministic across chunk sizes/thread
counts and exported to `<output>.offtarget-amplicons.tsv`. Only offending
candidate combinations are removed; an isolated strong primer hit does not
remove a pair. If no candidates remain for a target, screening fails explicitly
without relaxing the cutoff. `--off-target-mode single` preserves the prior
independent-primer filtering path. No reference still skips screening.

Scope: this mode checks retained intended F/R combinations, not all cross-target
primer combinations in the multiplex pool. The reference is treated as
background: intended loci are NOT automatically recognized or excluded.
An intended product present in -r is reported and rejected under the same
criteria; locus exclusion requires explicit mapping and is not implemented.
Seed-based hit discovery retains the original search's recall limitations.
Predicted compatible sites indicate potential amplification, not measured yield.
