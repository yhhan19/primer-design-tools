# Landscape and panel workflow

[Quick start](../README.md)

## Quick start: landscape → choose a cutoff → improve the panel → view SVG

### 1. Set input and create a fresh output directory

Run these commands from the repository root in Linux, macOS, or WSL:

```bash
targets=/path/to/aligned-targets
out=results/my-panel
mkdir -p "$out"
```

Replace `targets` with one aligned FASTA or a directory containing the alignments
you want to include. Every lowercase `.fa`/`.fasta` in a directory becomes a
separate target. Choose a fresh `out` to preserve previous reports. `-o` below
is a filename **prefix**, not a directory.

These examples omit `-r`, so they skip background screening. If using a background
reference, supply the same `-r` path in both runs and keep other design settings
consistent. Use your own aligned targets and include only the files you intend
to analyze in the input directory.

### 2. Create the landscape

```bash
./bin/dpro -i "$targets" -o "$out/landscape" \
  --landscape-only --dimer-landscape-points 10
```

This generates candidate primers and optimizes affected-primer count separately
at ten cutoffs per pool. It skips final panel selection and expansion. Defaults
sample evenly from −5000 to −2000 cal/mol; add `--dimer-landscape-adaptive` for
adaptive allocation. In an args file, set `dimer_landscape_points = 10`.

Open **`$out/landscape.landscape.svg`** in a browser. It shows both pools.
The numeric results are `$out/landscape.pool1.landscape.tsv` and
`$out/landscape.pool2.landscape.tsv`.

### 3. Inspect both pools and choose your cutoff

```bash
cutoff=-3000  # Example only: replace with the cutoff chosen according to the landscape.
```

The **min-max energy** is the worst interaction energy in the panel found by
the min-max search, which tries to make that worst interaction as high as
possible. At a cutoff equal to that panel's worst energy, it has zero affected
primers.

Optimizing only the single worst interaction can constrain the other primer
choices. Allowing a small number of affected primers can give the search more
flexibility to improve the rest of the panel and identify regions to revisit
through expansion.

A cutoff slightly above the best-found min-max energy can be a useful starting point for exploring expansion.

The horizontal axis is the energy cutoff; the vertical axis is affected-primer
percentage. An interaction counts as risky when its energy is **below**
the cutoff. Each participating primer counts once, even if it has several risky
partners.

The program does not choose a cutoff automatically. A more negative cutoff can
reduce the reported affected count simply by changing what counts as risky.
Choose according to the assay criteria, rather than selecting the smallest
plotted count. Search results are heuristic.

### 4. Rerun selection and expansion at that cutoff

```bash
./bin/dpro -i "$targets" -o "$out/final" \
  --dimer-dg-cutoff "$cutoff"
```

The defaults select `min-affected-primers` and enable one affected-PDR expansion
pass in the workflow. The program selects an initial panel,
expands affected regions, generates additional candidates and tries selection
again. A trial must pass layout checks and improve the dimer performance to be accepted; otherwise the original panel is retained.

Use the **same input and design settings** as the landscape run. If using args
files, use separate landscape and final presets.
Loading the saved `landscape.args` unchanged would retain landscape-only mode.
For a comparison without expansion, add `--no-expand-affected-pdrs` to a separate run.

Both output prefixes share `$out/.dpro-dimer-cache`, so unchanged interaction
energies are reused. Expansion calculates missing interactions for new candidates.
Upstream candidate-generation stages still run. If the runs use different output
directories, pass the same `--dimer-cache /path/to/shared-cache` to both.
Logs report reused and newly calculated interaction counts.

### 5. Check the outcome and view the SVGs

```bash
cat "$out/final.expansion.tsv"
cat "$out/final.before-expansion.dimer-summary.tsv"
cat "$out/final.dimer-summary.tsv"
```

Read `# status` in the expansion audit. `accepted` means the improved trial became
the final panel; other statuses retain the original panel. Compare before/after
counts at the same cutoff using these summaries. A landscape point may come from
a different search and is not the expansion baseline.

Open **`$out/final.coverage.svg`** in a browser to inspect the panel actually
selected. This SVG is generated automatically; no plotting command is needed.
Orange primer outlines identify below-cutoff interactions.

For a separate before/after view, generate the comparison when a trial exists:

```bash
if [ -f "$out/final.expanded.coverage.svg" ] && [ -f "$out/final.expanded.bed" ]; then
  python3 scripts/plot_pdr_expansion.py "$out/final" --no-links
fi
```

Open **`$out/final.expansion-comparison.svg`** if generated. It compares **Before**
with **Trial**, and its header states whether the trial was accepted. Gray dashed
boxes show original PDRs, purple spans show expanded regions, and red primer
marks identify affected primers. `--no-links` omits connection lines; omit that
flag to retain them. If no trial was generated, use `final.coverage.svg` and the
expansion audit instead.

From WSL, open the output folder and then open the SVGs in your browser:

```bash
explorer.exe "$(wslpath -w "$out")"
```

| File under `out` | What to inspect |
| --- | --- |
| `landscape.landscape.svg` | Both pools across sampled cutoffs; inspect before choosing a cutoff |
| `final.coverage.svg` | The actual final panel, whether expansion was accepted or not |
| `final.expansion-comparison.svg` | Before versus trial; generated separately when a trial exists |
| `final.expansion.tsv` | Expansion status and candidate-generation audit |
| `final.before-expansion.dimer-summary.tsv`, `final.dimer-summary.tsv` | Initial versus final counts at your chosen cutoff |
| `landscape.args`, `final.args` | Effective settings for reproducibility |
| `landscape.log`, `final.log` | Complete run logs |

See the sections below for landscape, expansion, cache and logging details,
and the [output reference](outputs.md) for report formats.

## Landscape details

`--dimer-landscape-points` defaults to 0 (off); enabled values are 2–100.
Each cutoff gets an independent affected-primer optimization. Landscape budgets
are `--dimer-landscape-restarts` (200) and `--dimer-landscape-iterations` (1000
proposals per restart). The default range is −5000 to −2000 cal/mol.

Uniform sampling is the default. Adaptive mode starts at the lower and upper
endpoints, then the configured selection cutoff if interior. It repeatedly
bisects the interval maximizing `width * (1 + endpoint affected-count change)`.
Discovered panels are reused as search seeds and improve best-known counts at
already sampled cutoffs.

The combined `.landscape.svg` shows both pools. Plotted counts are each cutoff's
own-search result (`raw_affected_primers`); the TSV `affected_primers` column can
include improvements from panels discovered at other cutoffs.

`--landscape-only` skips main selection and selected-panel reports. It requires
landscape points, full-pipeline input, expansion disabled and no main-panel seed
prefix. Expansion defaults off automatically in this mode; explicitly enabling
it is an error. Landscape-only writes the combined SVG and per-pool tables,
omitting per-pool SVGs and primer witnesses. Normal landscape runs also write
per-pool plots and `.landscape.primers.tsv` witnesses.

## Expansion and warm starts

Expansion defaults on for full midpoint panel runs. It defaults off for fixed
layout, PDR-only and landscape-only; explicit enablement there errors. Disable
it with `--no-expand-affected-pdrs` or `expand_affected_pdrs = 0` in args.

The single expansion pass grows affected windows subject to target boundaries
and the alternating-pool layout. Original windows remain contained; unaffected
windows stay fixed. New pairs undergo the configured generation, geometry,
alignment and optional background checks. Original candidates remain available.
Unchanged interaction energies are reused; the initial panel seeds reselection.

The trial is checked for permitted product lengths, window containment,
preservation of covered positions and same-pool layout. It is accepted only when
valid and improved under the configured objective without worsening either pool.
Otherwise the initial panel remains final. `.expanded.*` outputs describe the
trial even when rejected; `.coverage.svg` describes the actual final panel.
No-op expansion can produce no trial reports.

`--dimer-start-prefix` is an advanced warm start from
`<prefix>.poolN.<method>.primers.tsv`. Candidate indices, sequences, positions,
lengths and allowed combinations must match. This is not an interrupted-run
checkpoint and does not cache upstream stages. Changing the candidate set can
invalidate seeds. Persistent energy caching is separate and can still reuse
unchanged sequence interactions.

## Reusing interaction energies

`--dimer-cache auto` stores persistent energies in
`<output-directory>/.dpro-dimer-cache`. Different prefixes in the same directory
share it. Use an explicit directory across different output locations, or `off`
to disable reuse. The args key is `dimer_cache`.

Keys use ordered physical 5′→3′ primer sequences and temperature bounds. Cache
compatibility checks concentrations, actual loaded thermodynamic parameter
contents, and the compiled scoring implementation. Changing only a cutoff,
objective, seed, graph index or PDR boundary does not invalidate unchanged
sequence interactions. Nodes and allowed pair choices are rebuilt from current
candidates; cached energies are inserted into the rebuilt graph.

Each graph logs `Dimer cache: reused X, computed Y of Z interactions`. Complete
immutable batches are published atomically. Interrupted partial batches are
ignored; damaged published batches cause an error. Concurrent processes publish
separate batches safely, but can independently calculate interactions absent
from their initial snapshots. Cache files can be removed to force recalculation.

Existing graph exports and panel seed files are not imported. The first
cache-enabled run calculates missing energies. Makefile builds fingerprint the
scoring wrapper and linked libthal object; custom builds use the wrapper build
timestamp as implementation identity.

## Logs and console

Every normal run archives stdout/stderr to `<output-prefix>.log`, including
startup settings, progress, warnings and errors. Argument errors before a log
can be opened fall back to `dpro-error.log` in the working directory.

`--console-mode auto` uses a live view on a terminal and plain, append-only output
when redirected. Completed stages collapse to summaries; the complete history
remains in the plain-text log. `live` forces the terminal view; `plain` disables
it. Both modes produce the same design results. Use `--log-file /path/to/run.txt`
for a custom archive. Reusing a prefix or log path can overwrite previous files.

Full runs can be expensive: main selection defaults to 400 restarts per
objective/pool, with up to `5 * -I` proposals per restart. Landscape searches have
separate budgets. Graph exports and persistent caches can also be large.
