#!/usr/bin/env bash
# Run from WSL. Uses the repository's BTV alignments, without a background reference.
set -euo pipefail

btv_repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$btv_repo_root"

usage() {
  cat <<'EOF'
Usage:
  bash examples/btv_workflow.sh landscape [OUTPUT_DIR]
  bash examples/btv_workflow.sh design CUTOFF [OUTPUT_DIR]
  bash examples/btv_workflow.sh all CUTOFF [OUTPUT_DIR]

Default action: landscape
Default output directory: results/btv

Examples:
  bash examples/btv_workflow.sh landscape
  bash examples/btv_workflow.sh design -3000
  bash examples/btv_workflow.sh all -3000 results/btv-test

Inspect the landscape before choosing a cutoff. The 'all' action uses the
cutoff supplied on the command line without pausing or selecting it for you.
Use another output directory to preserve outputs from an earlier run.
EOF
}

btv_action="${1:-landscape}"
if (($#)); then shift; fi
btv_cutoff=''
case "$btv_action" in
  -h|--help|help) usage; exit 0 ;;
  landscape)
    if (($# > 1)); then usage >&2; exit 2; fi
    btv_output_dir="${1:-results/btv}"
    ;;
  design|all)
    if (($# < 1 || $# > 2)); then usage >&2; exit 2; fi
    btv_cutoff="$1"
    btv_output_dir="${2:-results/btv}"
    ;;
  *) usage >&2; exit 2 ;;
esac

if [[ ! -x bin/dpro ]]; then
  echo 'Build first: make -C primer3/src && make -j8' >&2
  exit 1
fi
if [[ ! -d scripts/aligned_sequences ]]; then
  echo 'Missing BTV input directory: scripts/aligned_sequences' >&2
  exit 1
fi
if [[ "$btv_action" != landscape ]]; then
  command -v python3 >/dev/null || { echo 'Python 3 is required for plotting.' >&2; exit 1; }
  python3 - "$btv_cutoff" <<'PY'
import math
import sys
try:
    cutoff = float(sys.argv[1])
    if not math.isfinite(cutoff):
        raise ValueError()
except ValueError:
    sys.exit('CUTOFF must be a finite number, for example -3000.')
PY
fi
mkdir -p -- "$btv_output_dir"

run_landscape() {
  ./bin/dpro --args examples/btv-landscape.args \
    -i scripts/aligned_sequences \
    -o "$btv_output_dir/landscape"
  printf '\nLandscape SVG: %s\nFull log: %s\n' \
    "$btv_output_dir/landscape.landscape.svg" "$btv_output_dir/landscape.log"
}

run_design() {
  ./bin/dpro --args examples/btv-expand.args \
    -i scripts/aligned_sequences \
    --dimer-dg-cutoff "$btv_cutoff" \
    -o "$btv_output_dir/design"

  # A no-op expansion has no reselection trial to compare.
  btv_expansion_status="$(awk -F '\t' '$1 == "# status" {print $2}' "$btv_output_dir/design.expansion.tsv")"
  if [[ "$btv_expansion_status" != no_new_* && \
        -f "$btv_output_dir/design.expanded.bed" && \
        -f "$btv_output_dir/design.expanded.coverage.svg" ]]; then
    python3 scripts/plot_pdr_expansion.py "$btv_output_dir/design" --no-links
    printf '\nComparison SVG: %s\n' "$btv_output_dir/design.expansion-comparison.svg"
  else
    printf '\nNo reselection trial was generated. Panel SVG: %s\n' "$btv_output_dir/design.coverage.svg"
  fi
  printf 'Expansion audit: %s\nBefore: %s\nAfter: %s\nFinal panel: %s\nFull log: %s\n' \
    "$btv_output_dir/design.expansion.tsv" \
    "$btv_output_dir/design.before-expansion.dimer-summary.tsv" \
    "$btv_output_dir/design.dimer-summary.tsv" \
    "$btv_output_dir/design.bed" "$btv_output_dir/design.log"
}

case "$btv_action" in
  landscape)
    run_landscape
    printf '\nInspect the landscape and choose a cutoff. Example next command:\n'
    printf 'bash examples/btv_workflow.sh design -3000 %q\n' "$btv_output_dir"
    ;;
  design) run_design ;;
  all) run_landscape; run_design ;;
esac
