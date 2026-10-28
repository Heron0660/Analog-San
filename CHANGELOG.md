# Changelog

## 1.7.0 - 2026-10-28

- Added optional deterministic MaxMin fingerprint selection within manual rule tiers for small direct-edit portfolios.
- Kept the existing tier precedence and tier counts; diversity selection is not described as a synthesis-feasibility predictor.
- Added a paired 28-prototype ablation at requested set sizes 5 and 10. The observed within-set nearest-neighbour similarity decreased, with a small trade-off in parent similarity at size 5 and a slight increase in SA Score; see `study/v1.7`.
- Added GUI control, export field, validation, and tests for the new selection mode.

## 1.6.0 - 2026-10-28

- Added a validated JSON catalogue for sourced transformation records, distinguishing reaction-class precedent from exact-transformation precedent.
- Candidate outputs and reports expose evidence status, source metadata, conditions, and route-verification status. Imported evidence never marks a substrate-specific route verified automatically.
- Clarified that manual edit tiers and composite scores are heuristics, not synthetic-feasibility predictions.
- Added tests for catalogue validation and evidence reporting.
- Re-ran fixed panels and the archived matched-input CReM comparison. Aggregate counts and comparison metrics were unchanged from v1.5.0; no comparative-performance improvement is claimed.
- The distributed evidence catalogue remains empty (0 source records across 18 registered operations).
