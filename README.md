# Analog-San

Analog-San is a desktop research tool for generating and prioritizing candidate molecular analogues from a user-supplied prototype. Its 18 registered graph edits include classical and non-classical transformations, including replacement of aromatic chlorine, bromine, or iodine by methyl, ethyl, isopropyl, or tert-butyl groups. It supports Lipinski and Veber filters, heuristic synthesis-priority ranking, molecular-property summaries, and CSV, SMILES, and MOL2 export.

## Scope and limitations

The software proposes structures; it does not predict binding affinity, biological efficacy, ADMET outcomes, toxicity, or experimentally validated synthetic routes. Priority and composite scores are heuristic ranking aids, not calibrated probabilities or evidence of synthesizability. The aromatic halide replacements are structure-enumeration rules, not reaction templates or route predictions. Every generated structure requires cheminformatics review and experimental validation. Filters can exclude useful compounds and should be selected deliberately.

## Installation and use

Requires Python 3.10 or newer and a working Tk installation. RDKit installation is platform-dependent; use the official RDKit installation instructions if `pip install` is unavailable for your platform.

```bash
python -m pip install -r requirements.txt
python main.py
```

To run the reproducible benchmark study, install `study/requirements-study.txt`, then follow `study/README.md`. The benchmark scripts use fixed random seeds and write their results under `study/dados/`.

## Source layout

- `main.py`, `gui_app.py`: application entry point and interface.
- `chem_engine.py`: analogue-generation and molecular-property logic.
- `validation_engine.py`: validation and figure-generation support.
- `study/`: benchmark inputs, analysis scripts, generated tables, figures, and run metadata.

## Reproducibility

The study scripts document their versions, seeds, and checksums in `study/metadata_run.json`. Re-run the scripts to regenerate benchmark outputs and publication figures. Results describe computational experiments on the included prototype set; they are not prospective validation.

## License and citation

Released under the MIT License; see `LICENSE`. Citation metadata is provided in `CITATION.cff`. Add the archival DOI and final author list after depositing a tagged release in a DOI-minting archive.

## Contact

Please use the repository issue tracker for reproducible bug reports. Do not submit confidential structures or unpublished patent material in public issues.
