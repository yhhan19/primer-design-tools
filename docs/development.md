# Validation and repository layout

[Quick start](../README.md)

## Validation and repository layout

After building, run the regression suites from Linux/WSL:

```bash
python3 tests/run_regressions.py
python3 tests/run_regressions.py --cpp --jobs 2
```

See [testing-workflow.md](testing-workflow.md) for test scope, logs and
reproducibility checks. Cache unit tests cover concurrent requests, persistence,
invalidation and corruption checks; the integration test covers landscape/final
reuse, expansion and uncached equivalence. Test availability does not imply a
particular checkout has passed every test.

- `src/`: CLI, pipeline, PDR layout, pair scoring, screening, graph selection and cache.
- `primer3/`: Primer3 submodule and thermodynamic parameters.
- `examples/`: saved presets and workflow helpers.
- `docs/`: workflow and testing guides.
- `tests/`: CLI and C++ regression checks.
- `scripts/`: plotting, analysis and example datasets; some scripts need extra dependencies.
- `data/`: example inputs and historical outputs.
- `build/`, `bin/`: build products.

Coverage and interaction reports describe computational criteria. Alignment,
background screening and report generation remain separate from experimental
validation.
