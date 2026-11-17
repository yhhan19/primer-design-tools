# WSL testing workflow

Run these commands in WSL, from this checkout. The real-data examples omit `-r`
and therefore do not run background off-target screening.

## 1. Build and run the regression suite

```bash
cd /mnt/c/work/off-target-binding/primer-design-tools
make -C primer3/src
make -j8

# CLI and Python regression tests:
python3 tests/run_regressions.py

# Include all C++ regression programs:
python3 tests/run_regressions.py --cpp --jobs 2
```

Results and per-test logs are retained in `build/software-tests/`; the summary is
`build/software-tests/summary.json`. The runner returns nonzero if any test fails.
Build again before testing if you changed C++ source or headers. C++ tests use
assertions and should not be compiled with `-DNDEBUG`.

To rerun a targeted group:

```bash
python3 tests/run_regressions.py --cpp --only \
  test_candidate_pairs test_pair_msa test_expansion_layout test_run_logging
```

The suite covers:

- CLI validation, saved args, defaults, and disabled/removed options.
- Primer orientation, dimer temperature ranges and concurrent calculations.
- Pair coordinates/product lengths, allowed-pair selection and MSA retention.
- Joint sequence coverage, 3-prime mismatches, ambiguous sites and ranking ties.
- Expansion clipping, unaffected PDRs remaining fixed, same-pool separation,
  acceptance/fallback, and final report consistency.
- Landscape-only mode, cutoff decision/result logs, combined SVG contents and
  absence of selected-panel outputs.
- Full text logs, live terminal collapse, plain fallback, errors and interruption.
- Small synthetic off-target controls in the C++ suite. These do not load a
  livestock or human reference genome.
- Cache unit tests for concurrent requests, disk reload, invalidation, concurrent
  batch publication, and corruption checks. Graph reuse integration checks cover
  fresh-process landscape-to-selection reuse, expansion, and uncached parity.

### Interaction cache reuse

Use the same output directory for landscape and final runs to share the default
`.dpro-dimer-cache`, or supply the same `--dimer-cache` directory to both runs.
Changing the cutoff does not invalidate energies. Expansion reuses old
interactions and calculates missing ones. `--dimer-cache off` provides an
uncached comparison; compare final panels and weighted edges, not timings/logs.

## 2. Small real-data acceptance test: FMD

This input has 28 sequences aligned over 8,180 positions. It is a useful quick
check of the complete selection/expansion pipeline, without a background scan.

```bash
mkdir -p results/validation

./bin/dpro --args examples/btv-expand.args \
  -i scripts/FMD.msa.fasta.filtered.fasta \
  --dimer-dg-cutoff -3000 \
  -o results/validation/fmd

cat results/validation/fmd.expansion.tsv
cat results/validation/fmd.before-expansion.dimer-summary.tsv
cat results/validation/fmd.dimer-summary.tsv

# Run this when a reselection trial was generated:
python3 scripts/plot_pdr_expansion.py results/validation/fmd
```

The BTV-named args file is a design preset; `-i` overrides its input. The final
effective settings are saved as `results/validation/fmd.args`. A previous run
using 400 restarts and `-I 1000` retained 25 amplicons and improved 2 affected
primers to 0. Treat that as a regression reference for this input, version and
seed, not a general quality guarantee.

Open `results/validation/fmd.expansion-comparison.svg`. Gray boxes are original
PDRs, purple spans are expansions, and red arrows/underlines are affected primers.

## 3. Larger alignment smoke test: Zika

This input has 372 sequences aligned over 10,753 positions. Use a smaller main
search budget for a functional smoke test:

```bash
./bin/dpro --args examples/btv-expand.args \
  -i scripts/ZIKA.msa.fasta.filtered.fasta \
  --dimer-dg-cutoff -3000 \
  --dimer-restarts 200 -I 200 \
  -o results/validation/zika-smoke

cat results/validation/zika-smoke.expansion.tsv
cat results/validation/zika-smoke.before-expansion.dimer-summary.tsv
cat results/validation/zika-smoke.dimer-summary.tsv
python3 scripts/plot_pdr_expansion.py results/validation/zika-smoke
```

`-I 200` means 1,000 proposals per restart; 200 restarts is 200,000 proposals
per pool/search. This is a smoke budget, not an optimality guarantee. The tested
run selected 35 amplicons and accepted an improvement from 2 affected primers to
0, while preserving coverage and same-pool separation.

For the actual cutoff-choice workflow, first run the landscape separately:

```bash
./bin/dpro --args examples/btv-landscape.args \
  -i scripts/ZIKA.msa.fasta.filtered.fasta \
  --dimer-landscape-points 10 \
  --dimer-landscape-restarts 100 \
  --dimer-landscape-iterations 500 \
  -o results/validation/zika-landscape
```

Open `results/validation/zika-landscape.landscape.svg`. This step produces no
selected primer panel. Choose a cutoff, then use the expansion command above
with that cutoff. A lower-budget landscape may report more affected primers
than a stronger final design search.

## 4. Reproducibility: same seed and settings

After the FMD run:

```bash
./bin/dpro --args results/validation/fmd.args \
  -S 42 -o results/validation/fmd-repeat

cmp results/validation/fmd.bed results/validation/fmd-repeat.bed
```

`cmp` should print nothing and exit 0 when the panels are identical. Compare BED
files, not logs or args timestamps. Change the seed for independent exploration;
a different seed is allowed to produce a different valid panel.

## 5. No-affected-primer behavior

Use an extreme cutoff only as a software edge case:

```bash
./bin/dpro --args examples/btv-expand.args \
  -i scripts/FMD.msa.fasta.filtered.fasta \
  --dimer-dg-cutoff -1000000 \
  --dimer-restarts 20 -I 20 \
  -o results/validation/fmd-no-affected

cat results/validation/fmd-no-affected.expansion.tsv
cmp results/validation/fmd-no-affected.before-expansion.bed \
    results/validation/fmd-no-affected.bed
```

Expect zero affected primers, `no_new_compatible_pairs`, and an unchanged panel.
Do not interpret this numerical test cutoff as a suitable design threshold.
There is no reselection trial to plot when expansion is a no-op.

## 6. Expansion-radius sensitivity

```bash
for radius in 5 20 40; do
  ./bin/dpro --args examples/btv-expand.args \
    -i scripts/FMD.msa.fasta.filtered.fasta \
    --dimer-dg-cutoff -3000 --pdr-expansion-bases "$radius" \
    -o "results/validation/fmd-expand-${radius}"
done
```

Compare each run's initial and final dimer summaries at the same cutoff. Check
that expanded windows contain the originals and remain within the specified
radius and layout constraints. A larger radius is not guaranteed to give a
better heuristic result. Rejected or tied revisions must leave the original
panel in the final `.bed`.

## 7. Scale test: COVID subset

This input has 342 sequences aligned over 29,903 positions. Run it last because
the candidate graph can dominate runtime and memory even with a small SA budget.

```bash
time ./bin/dpro --args examples/btv-landscape.args \
  -i scripts/covid.aligned.fasta.zst.subset.fasta \
  --dimer-landscape-points 6 \
  --dimer-landscape-restarts 50 \
  --dimer-landscape-iterations 300 \
  -o results/validation/covid-landscape
```

Despite its name, this input is an ordinary FASTA file. This command is a scale
smoke test, not a recommendation to use its search budget for final design.
Do not compare affected fractions across different datasets as if they measured
software quality; candidate availability and sequence variation differ.

## What to check in each design run

1. The input, seed, cutoff, length requirements and search budgets in `.args`.
2. The expansion status before interpreting any trial as the final panel.
3. Affected primers before/after at the same cutoff, plus newly affected primers.
4. Actual product lengths, zero same-pool overlap, and no loss of original
   covered positions. The regression tests check these independently.
5. MSA scores for selected replacements in `.pair-msa.tsv`: MSA filtering is
   relative and does not guarantee that the replacement has non-decreasing MSA
   robustness. Geometric coverage is different from MSA sequence coverage.
6. Warnings/errors in `<prefix>.log`; live console collapse does not remove them
   from that file. Use a new prefix to preserve earlier outputs.

Computational constraints and these regression results do not establish
experimental amplification performance.
