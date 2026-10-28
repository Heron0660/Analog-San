# Reproducible benchmark study

This folder contains the inputs and scripts supporting the computational study. It is a reproducibility package, not an independent validation of pharmacological activity or experimental synthesis.

## Environment

Use Python 3.13.9 with the packages pinned in `requirements-study.txt`. RDKit builds can be platform-specific; install the matching official build if the pinned package is unavailable.

## Re-run

From this directory:

```bash
python scripts/run_benchmarks.py
python scripts/analyze.py
python scripts/make_figures.py
```

The benchmark script regenerates the data tables and run metadata. The analysis script regenerates summary outputs. The figure script writes publication figures in PNG and SVG formats to `figures/`. Fixed seeds are recorded in `metadata_run.json`.

The uniform-sampling comparator changes only the random exploration choice weights. Direct-stage ordering and priority-tier definitions remain identical to the weighted method. Because those tier definitions also contribute to ranking, this is an algorithmic sampling comparison, not a blinded comparison against an independent synthesis planner.

All scores, alerts, and synthesis-priority labels are heuristic. No binding, efficacy, ADMET, toxicity, or route feasibility predictions are validated here.
