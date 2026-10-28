"""Reanalyse the archived CReM candidate pools against Analog-San 1.4."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger, rdBase

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import chem_engine as ce
import validation_engine as ve
from run_crem_comparison import PROTOTYPES, SEEDS, TARGET

RDLogger.DisableLog("rdApp.*")
SA_LIMIT = 6.0
REACTIVE_TOKENS = ("reativo", "epóxido tensionado", "energético", "michael acceptor")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def common_filter(mol: Chem.Mol) -> tuple[bool, dict]:
    props = ce.calculate_properties(mol)
    sa = ve.synthetic_accessibility_score(mol)
    alerts = ve.chemical_alerts(mol)
    reactive = [
        alert for alert in alerts
        if any(token in alert.casefold() for token in REACTIVE_TOKENS)
    ]
    accepted = (
        ce.lipinski_filter(props)
        and ce.veber_filter(props)
        and sa <= SA_LIMIT
        and not reactive
    )
    return accepted, {**props, "SAScore": sa, "AlertCount": len(alerts),
                      "AlertNames": "; ".join(alerts), "ReactiveAlerts": "; ".join(reactive)}


def describe(name: str, seed: int, tool: str, smi: str) -> dict:
    mol = Chem.MolFromSmiles(smi)
    parent = Chem.MolFromSmiles(PROTOTYPES[name])
    accepted, descriptors = common_filter(mol)
    assert accepted
    return {
        "Prototype": name, "Seed": seed, "Tool": tool,
        "SMILES": Chem.MolToSmiles(mol, isomericSmiles=True),
        **descriptors,
        "SimilarityTanimoto": ve.tanimoto_similarity(parent, mol),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    pd.DataFrame(rows).to_csv(path, index=False, encoding="utf-8-sig")


def run(old_crem_pool: Path, old_run_metadata: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    old_pool = pd.read_csv(old_crem_pool)
    previous_runs = pd.read_csv(old_run_metadata)
    crem_records: list[dict] = []
    analog_records: list[dict] = []
    run_rows: list[dict] = []

    for name, parent_smiles in PROTOTYPES.items():
        for seed in SEEDS:
            old_group = old_pool[(old_pool.Prototype == name) & (old_pool.Seed == seed)]
            eligible = []
            for smi in old_group.SMILES:
                mol = Chem.MolFromSmiles(smi)
                if mol is not None and common_filter(mol)[0]:
                    eligible.append(Chem.MolToSmiles(mol, isomericSmiles=True))
            eligible = sorted(set(eligible))
            selected = set(random.Random(seed).sample(eligible, min(TARGET, len(eligible))))
            for smi in eligible:
                row = describe(name, seed, "CReM", smi)
                row["SelectedForTarget"] = int(smi in selected)
                crem_records.append(row)

            generated = ce.generate_analogs(
                parent_smiles,
                TARGET,
                use_veber=True,
                selected_routes=[ce.ROUTE_CLASSIC, ce.ROUTE_NON_CLASSIC],
                random_seed=seed,
                max_sa_score=SA_LIMIT,
                allow_multistep=False,
            )
            for item in generated:
                analog_records.append(describe(name, seed, "Analog-San", item["SMILES"]))

            old_meta = previous_runs[
                (previous_runs.Prototype == name)
                & (previous_runs.Seed == seed)
                & (previous_runs.Tool == "CReM")
            ].iloc[0]
            run_rows.extend([
                {"Prototype": name, "Seed": seed, "Tool": "CReM",
                 "PoolBeforeCommonFilters": len(old_group),
                 "PassCommonFilters": len(eligible), "Returned": len(selected),
                 "Target": TARGET, "OldRawMutate": old_meta.RawMutate,
                 "OldRawGrow": old_meta.RawGrow,
                 "Source": "re-filtered archived full candidate pool"},
                {"Prototype": name, "Seed": seed, "Tool": "Analog-San",
                 "PoolBeforeCommonFilters": "NA", "PassCommonFilters": len(generated),
                 "Returned": len(generated), "Target": TARGET,
                 "OldRawMutate": "NA", "OldRawGrow": "NA",
                 "Source": f"Analog-San {ve.APP_VERSION} direct-only generator; 18 registered operations"},
            ])

    write_csv(output_dir / "records_crem_comparison.csv", crem_records)
    write_csv(output_dir / "records_analogsan_comparison.csv", analog_records)
    write_csv(output_dir / "run_metadata_comparison.csv", run_rows)

    metadata = {
        "protocol": "Matched-input reanalysis: 28 prototypes, 3 fixed seeds, target 100.",
        "software": f"Analog-San {ve.APP_VERSION} (18 registered operations)", "validation": ve.VALIDATION_VERSION,
        "evidence_catalog": "empty; no reaction records catalogued; no route verification claimed",
        "RDKit": rdBase.rdkitVersion, "Python": sys.version.split()[0],
        "CReM": "0.3.2 candidate pool generated previously with archived protocol",
        "common_filters": {
            "Lipinski": "MW<=500, MolLogP<=5, HBD<=5, HBA<=10; zero violations",
            "Veber": "rotatable bonds<=10 and TPSA<=140 A^2",
            "SA_score_max": SA_LIMIT,
            "reactive_alert_exclusions": list(REACTIVE_TOKENS),
            "PAINS_Brenk_other_alerts": "reported, not hard-excluded",
        },
            "Analog-San": "Current 18-operation direct prototype registry; multi-step exploration disabled; same common filters; up to 100 outputs. This is a post hoc descriptive reanalysis, not an independent validation of the newly added transformations.",
        "CReM": "Full archived pools re-filtered with the common filters, then seeded uniform sample without replacement capped at 100.",
        "CReM_database_bias": "The archived CReM candidate pools derive from a ChEMBL33 fragment database prefiltered for source SA<=2 and structural alerts; this remains a comparator limitation.",
        "seeds": list(SEEDS),
        "inputs": {
            old_crem_pool.name: sha256(old_crem_pool),
            old_run_metadata.name: sha256(old_run_metadata),
        },
        "outputs_sha256": {},
    }
    for filename in ("records_crem_comparison.csv", "records_analogsan_comparison.csv", "run_metadata_comparison.csv"):
        metadata["outputs_sha256"][filename] = sha256(output_dir / filename)
    (output_dir / "metadata_comparison.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--crem-pool", type=Path, default=ROOT / "study" / "dados" / "records_crem_comparison.csv")
    parser.add_argument("--run-metadata", type=Path, default=ROOT / "study" / "dados" / "run_metadata_comparison.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "study" / "dados_v1.4")
    args = parser.parse_args()
    run(args.crem_pool.resolve(), args.run_metadata.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()
