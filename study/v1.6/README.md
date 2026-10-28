# Analog-San 1.6.0 reproducibility outputs

These files reproduce the fixed-panel characterization and matched-input descriptive comparison reported for Analog-San 1.6.0. The direct generator is deterministic under the evaluated settings; repeated seeds check invariance and are not independent samples.

## Reproduce

Install the study and comparison dependencies listed in `../requirements-study.txt` and `../requirements-comparison.txt`. From the repository root, run:

```bash
python study/v1.6/scripts/rerun_panels.py --output-dir study/v1.6/outputs/main
python study/v1.6/scripts/analyze_panels.py --data-dir study/v1.6/outputs/main
python study/v1.6/scripts/rerun_comparison.py --crem-pool study/dados/records_crem_comparison.csv --run-metadata study/dados/run_metadata_comparison.csv --output-dir study/v1.6/outputs/comparison
python study/v1.6/scripts/summarize_comparison.py --data-dir study/v1.6/outputs/comparison --figure-dir study/v1.6/outputs/figures
python study/v1.6/scripts/make_figures.py --data-dir study/v1.6/outputs/main --figure-dir study/v1.6/outputs/figures
python study/v1.6/scripts/catalog_census.py --output study/v1.6/outputs/catalog_census.csv
```

The comparison re-filters archived CReM candidate pools; it does not regenerate CReM outputs. The pool derives from ChEMBL33 fragments prefiltered for source SA≤2 and structural alerts, confounding SA and alert contrasts. Results are descriptive and do not assess biological activity, novelty, expert-rated quality, or experimental synthesis feasibility. The temporal challenge in the archived v1.5 analysis is outcome-informed and post hoc; it is not independent validation.

## Evidence catalogue

The shipped `evidence_catalog.json` is empty. The generated census therefore reports 0 catalogued source records and 0 route-verified transformations across 18 registered operations. This describes curation status only; it does not imply that literature precedents are absent. See [`../../docs/evidence_catalog.txt`](../../docs/evidence_catalog.txt).
