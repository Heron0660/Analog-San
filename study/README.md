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

## CReM comparative benchmark

The paired comparison uses CReM 0.3.2 and the ChEMBL33 SA<=2 fragment database archived at Zenodo DOI `10.5281/zenodo.16909329`. Download `chembl33_sa2_f5.db.gz` from the record, verify its SHA-256 against `metadata_comparison.json`, and extract the SQLite file. The database is filtered for source-molecule SA<=2 and structural alerts, so SA and alert comparisons are descriptive and confounded by database construction.

Use Python 3.12 and `requirements-comparison.txt` for this comparison. From `study/`:

```bash
python scripts/run_crem_comparison.py --db path/to/chembl33_sa2_f5.db
python scripts/analyze_crem_comparison.py
```

The comparison uses 28 prototypes, three fixed seeds, the same strict Lipinski and Veber filters, CReM `mutate`/`grow` settings, a cap of 1,000 candidates per CReM mode and seed, and a seeded uniform sample of up to 100 passing structures per run. The target was increased from 30 because the lower target saturated many runs and obscured output capacity. The complete filtered CReM candidate pools, selected records, per-run metadata, and derived summaries are provided in `dados/`. This is an exploratory matched-input software comparison, not an independent synthesis-feasibility or biological validation.
