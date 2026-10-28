# Changelog

## 1.6.0 - 2026-10-28

- Added a validated JSON catalogue for sourced transformation records, distinguishing reaction-class precedent from exact-transformation precedent.
- Candidate outputs and reports expose evidence status, source metadata, conditions, and route-verification status. Imported evidence never marks a substrate-specific route verified automatically.
- Clarified that manual edit tiers and composite scores are heuristics, not synthetic-feasibility predictions.
- Added tests for catalogue validation and evidence reporting.
- Re-ran fixed panels and the archived matched-input CReM comparison. Aggregate counts and comparison metrics were unchanged from v1.5.0; no comparative-performance improvement is claimed.
- The distributed evidence catalogue remains empty (0 source records across 18 registered operations).
