# CLI reference and pipeline

[Quick start](../README.md)

## Inputs and pipeline

- **Targets (`-i`, required):** one FASTA alignment or a nonrecursive directory
  of alignments. Paths are sorted before assigning segment indices. Sequences
  within each alignment must have equal lengths; different segments may differ.
  A single sequence can also serve as an alignment.
- **Background (`-r`, optional):** reference FASTA for off-target screening.
  Intended target products are not automatically excluded from it.
- **Output (`-o`, required):** filename prefix; create its parent directory first.
- **Precomputed risk (`-a`, advanced):** one score per aligned position. This
  selects PDR-only mode. With directory input, the same file applies to each
  segment and must match every alignment's length.

FASTA input accepts DNA, IUPAC ambiguity codes and `-`, normalizes lowercase,
and requires a header before sequence data.

The full pipeline:

1. Build a consensus template and positional variation risk for each alignment.
2. Choose PDR placement (`midpoint` by default; `fixed` remains advanced).
3. Generate Primer3 candidates, check pair geometry and product length, then
   score pairs against the alignment and retain the highest-ranked pairs.
4. Optionally screen retained candidates against the background reference.
5. Alternate successive amplicons per segment between two pools and optimize
   interactions while preserving validated pair combinations.
6. By default, expand affected PDRs once and rerun selection, retaining the initial
   panel unless the trial passes validation and improves the configured objective.
7. Write the selected panel, coverage reports and dimer annotations.

## Main options and defaults

Use `--help` for common options and `--help-all` for the complete inventory.
Energies are in **cal/mol**, temperatures in **°C**, lengths in **bp**.

| Area | Options | Default |
| --- | --- | --- |
| Product size | `--amplicon-min` (`-Lx`), `--amplicon-max` (`-Ln`) | 252–420 |
| Layout | `--pdr-length` (`-Lp`), `--pdr-search` | 40, `midpoint` |
| Candidate feasibility | `--pdr-min-candidates` | 0 (adaptive) |
| PDR score bounds | `-Un`, `-Ux` | 0.1, 10000 |
| Primer size | `--p3-min-size`, `--p3-opt-size`, `--p3-max-size` | 18, 20, 25 |
| Primer Tm | `--p3-min-tm`, `--p3-opt-tm`, `--p3-max-tm` | 57, 60, 63 |
| Primer GC | `--p3-min-gc`, `--p3-max-gc` | 40–60% |
| Shortage warning / single-mode retention floor | `--candidate-min-return` | 5 |
| Pair retention | `--pair-keep-percent`, `--pair-min-return` | 20%, at least 5 when available; boundary ties retained |
| Pair coverage scoring | `--pair-max-mismatches`, `--pair-max-3prime-mismatches`, `--pair-3prime-window` | 2 total, 0 terminal, 5-base window |
| Off-target mode | `--off-target-mode` | `paired` |
| Off-target product range | `--off-target-amp-min`, `--off-target-amp-max` | 1–1000 |
| Search seed / mismatch limit | `-k`, `-H` | 8, 2 |
| Off-target energy cutoff | `-G` | −15000 |
| Workers | `-t` | 8 |
| Reference chunk / I/O block | `-C`, `-B` | 100000000, 134217728 |
| Objective | `--dimer-methods` | `min-affected-primers` |
| Dimer cutoff | `--dimer-dg-cutoff` | −3000 for cutoff objectives |
| Dimer temperature | `--dimer-temp-min`, `--dimer-temp-max` | 58–62 |
| Main search | `--dimer-restarts`, `-I`, `-S` | 400, 1000, 42 |
| Expansion | `--expand-affected-pdrs`, `--pdr-expansion-bases` | On for full midpoint runs, 20 bases per edge |
| Cache | `--dimer-cache` | `auto` |
| Concentrations | `--mv`, `--dv`, `--dntp`, `--dna` | 50, 1.5, 0.6 mM; 200 nM DNA |
| Non-dimer temperature | `--temp` | 37 |

The historical short names are unusual: `-Lx` is the **minimum**, `-Ln` the
**maximum**. Prefer their long aliases. `-H2` opts into an adjacent-mismatch
metric rather than standard `-H`. Pass signed energies, for example `-G -15000`.

## Pair ranking and retention

The mismatch settings assess the individual aligned sequences, not mismatches
against the consensus used to generate primers. Both primers must satisfy the
total and terminal mismatch limits on the same aligned sequence for it to count
as jointly covered. Gaps and ambiguous binding sites count as unknown rather
than matches; coverage fractions use all alignment sequences.

Pairs rank by joint coverage, then joint exact terminal-match fraction, then
weighted mismatch. Retention keeps the requested top percentage, at least the
minimum number when available, including all boundary ties. Individual primers
unused by retained pairs are removed. There is no absolute minimum-coverage
rejection threshold. Legacy individual-primer variation scoring is removed.

Primer3 constraints remain active. Repeated Primer3 validation of complete
pairs is disabled; pair geometry, actual product-length checks in midpoint mode,
and alignment ranking apply after candidate generation.

## Dimer objectives

Higher interaction energy is better. Objectives exclude self-dimers.

- `min-affected-primers`: minimize affected-primer count only; default, no
  secondary tie-breakers.
- `min-affected-amplicons`: minimize affected targets, then total cutoff deficit.
- `min-below-cutoff`: minimize below-cutoff interactions.
- `sum`, `min-max`, `top-k-bound`: advanced energy comparisons. They can be
  requested together. `top-k-bound` uses `--dimer-bad-percent` (default 1%) to set
  the allowed number of worst interactions; it is not an affected-primer percentage.

The three cutoff objectives must each be requested alone. An explicit
`--dimer-dg-cutoff` sets only the threshold; incompatible objectives cause an
error instead of an implicit method switch. For multi-method comparisons,
`--dimer-select` chooses the final panel; otherwise the first listed method wins.
Each method resets to `-S`. Equal dimer temperature bounds select one temperature;
`--temp` controls non-dimer scoring independently.

`--dimer-rounds` remains reserved: only 1 is accepted. Main selection uses up to
`5 * -I` proposals per restart. To reproduce the former main-search budget,
explicitly request `--dimer-restarts 4000`.

## Saved configurations and compatibility

Each parsed run writes effective settings to `<prefix>.args`. Files contain
`key = value` lines with internal names such as `len_amp`, `len_amp_min` and
`num_return`. Values use the rest of the line, supporting paths with spaces;
`#` begins a comment. Settings apply in command-line order; later values win,
including dimer methods, selection and cutoff.

```bash
./bin/dpro --args results/my-panel/final.args \
  -o results/my-panel/repeat -S 123
```

A separately supplied selection or cutoff stays in effect until overridden;
changing methods does not silently discard incompatible settings. Omit
`dimer_dg_cutoff` from configurations for non-cutoff comparison objectives.

`--candidate-min-return` replaces `--p3-num-return`; the latter remains a
deprecated alias with a warning. Existing `num_return` config keys still load.
It does not control Primer3 enumeration.

`--p3-weighted-keep-percent` is removed. Its legacy config key,
`p3_weighted_keep_percent`, warns and is ignored; new saved files omit it.
Other unknown options and keys remain errors, including removed panel-repair
keys. See the [decision record](../OPTION_DECISIONS.md) and
[option review](../OPTION_REVIEW.md).
