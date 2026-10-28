# Analog-San 1.7.1: small-set portfolio selection ablation

This study tests whether the optional tier-constrained MaxMin selector reduces redundancy in small direct-edit candidate sets relative to the default tier-priority ordering. It is an internal ablation, not a comparison of synthesis feasibility, activity, or general generator superiority.

## Reproduce

From the repository root, install `requirements-study.txt` and `requirements-comparison.txt`, then run:

```bash
python study/v1.7/scripts/benchmark_portfolio_selection.py --output-dir study/v1.7/outputs
```

The analysis uses the same 28 drug-like prototypes as the archived CReM comparison, with common Lipinski/Veber/SA/reactive-alert filters. It compares deterministic default and tier-diverse selections at requested sizes 5 and 10. Bootstrap intervals resample prototypes as paired units (20,000 replicates); structures or seeds are not treated as independent samples.

Within-set redundancy is the mean, over candidates, of each candidate's maximum Morgan radius-2, 2048-bit Tanimoto similarity to another candidate. Lower values mean lower nearest-neighbour similarity in this selected set. This metric does not establish scaffold novelty or pharmacological diversity. Manual priority-tier labels are retained and are not validated measures of synthetic ease.

## Results and limitations

In the archived evaluation, the selector reduced this redundancy metric by a mean 0.02782 at size 5 (95% paired prototype-bootstrap interval -0.03934 to -0.01662) and 0.00834 at size 10 (-0.01433 to -0.00295). At size 5, mean parent similarity decreased by 0.01285 and mean SA Score increased by 0.00950. At size 10, mean parent similarity decreased by 0.00151 (interval includes zero) and mean SA Score increased by 0.00521. Tier-count distributions did not change. These are small-set selection trade-offs, not evidence of synthetic feasibility or superiority over CReM.

The historical results are included in `outputs/`. Re-running with a different RDKit release may produce small descriptor or fingerprint differences; metadata records software versions and input hashes.
