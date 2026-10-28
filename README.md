# Analog-San

Analog-San 1.6.0 is an open-source desktop application for rule-based molecular analogue generation and computational triage. It applies 18 explicit graph edits, property filters, molecular descriptors, structural-alert flags, and reproducible structure/data exports.

## Scope

The software proposes structures; it does not predict binding, efficacy, ADMET, toxicity, novelty, or experimental synthesis success. Rule-priority tiers, SA Score, and composite scores are uncalibrated heuristics, not probabilities or route-feasibility estimates. Candidate structures may be impractical or unstable and require independent chemical review.

## Transformation evidence catalogue

`evidence_catalog.json` stores curated source records for registered transformations. Records are classified as `reaction_class` or `exact_transform`; neither classification verifies a route for a user's starting material. Exact-transformation records require reaction SMILES and an identified reviewer. Route verification remains false unless a complete substrate-specific route has been independently assessed. The released catalogue is empty: no records are currently curated for the 18 transformations.

Add records only after checking primary sources, structural correspondence, limitations, and reuse terms. See [`docs/evidence_catalog.txt`](docs/evidence_catalog.txt). The software validates the JSON structure and exposes evidence status and source metadata in candidate exports and reports; it does not infer evidence from a similar reaction.

## Installation

Requires Python 3.10 or newer, Tk, and RDKit. RDKit installation can be platform-specific.

```bash
python -m pip install -r requirements.txt
python main.py
```

Run the automated checks with:

```bash
python -m unittest discover -s tests -v
python main.py --smoke-test
```

## Reproducibility study

Benchmark inputs, scripts, fixed seeds, software versions, and generated outputs are provided under `study/`. The study is descriptive: it does not establish activity, synthetic feasibility, or comparative superiority. See [`study/README.md`](study/README.md) for protocols and caveats.

## Citation and license

The source is released under the MIT License. Citation metadata is in [`CITATION.cff`](CITATION.cff). This GitHub release has no DOI; an archival DOI may be added after deposition.
