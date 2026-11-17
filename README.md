# Primer Design Tools (`dpro`)

Tiled amplicon sequencing uses overlapping PCR products to cover a target.
Designing a multiplex panel requires accounting for sequence variation at primer
binding sites and interactions between primers sharing a reaction. `dpro`
combines variant-aware primer selection with thermodynamic interaction scoring,
then helps users explore dimer cutoffs, revisit primer design regions, and inspect the resulting panel in SVG reports.

BioRxiv preprint: Han and Boucher,
[*Scaling Variant-Aware Multiplex Primer Design*](https://doi.org/10.64898/2026.02.03.703607).

## Quick start

Run from the repository root in Linux, macOS, or WSL. You need Git, GNU Make
and a C++17 compiler; Python 3 is needed only for optional plots and tests.

### 1. Build and set your input

```bash
git submodule update --init --recursive
make -C primer3/src
make -j2

targets=/path/to/aligned-targets
out=results/my-panel
mkdir -p "$out"
```

Replace `targets` with one aligned FASTA or a directory containing only intended alignments. Choose a fresh output directory for each workflow.

### 2. Create a landscape (optional) and choose a cutoff

```bash
./bin/dpro -i "$targets" -o "$out/landscape" \
  --landscape-only --dimer-landscape-points 10
```

Open **`$out/landscape.landscape.svg`** in a browser and inspect both pools.
The **cutoff** is the predicted primer–primer interaction energy (ΔG, in cal/mol)
below which an interaction is labeled risky. More negative energies mean
stronger predicted binding.
A primer is **affected** if it has at least one risky interaction in its pool;
the landscape shows the percentage of affected primers at each sampled cutoff.

Choose a cutoff according to the assay criteria; the program does not select
it automatically. For exploring improvement, look near the transition from zero
to a modest affected fraction, before most primers become affected. A cutoff
slightly **above** (less negative than) the best-found min-max energy is a useful
starting point: it leaves risky interactions for expansion to address. This is
a search heuristic, not a validated assay threshold. See
[cutoff selection guidance](docs/workflow.md#choosing-an-informative-cutoff).

### 3. Run selection and try to improve the panel

```bash
cutoff=-3000  # Example: replace with your chosen cutoff.
./bin/dpro -i "$targets" -o "$out/final" --dimer-dg-cutoff "$cutoff"
```

This uses affected-primer selection and one expansion pass by default. A valid,
improved trial becomes final; otherwise the original panel is retained.
Use the same inputs and design settings as the landscape run. These examples
skip background screening; to enable it, add the same `-r reference.fasta` to both
commands. Unchanged interaction energies are reused automatically.

### 4. View the result

Open **`$out/final.coverage.svg`** for the actual final panel. Check
**`$out/final.expansion.tsv`** for whether expansion was accepted.

For a Before-versus-Trial SVG, follow the
[comparison plotting instructions](docs/workflow.md#5-check-the-outcome-and-view-the-svgs).
Full logs and effective settings are saved as `<prefix>.log` and `<prefix>.args`.

## Details

- [Installation and build](docs/installation.md)
- [Workflow, cutoff interpretation, SVG comparison and cache reuse](docs/workflow.md)
- [CLI options, defaults, pair scoring and configuration](docs/cli-reference.md)
- [Output files and formats](docs/outputs.md)
- [Testing](docs/testing-workflow.md) and [development](docs/development.md)

Use `./bin/dpro --help` for common options or `--help-all` for the complete list.
