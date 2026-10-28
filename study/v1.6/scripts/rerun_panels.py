"""Run the v1.4 direct-only main and extended prototype panels."""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger, rdBase

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import chem_engine as ce
import validation_engine as ve

RDLogger.DisableLog("rdApp.*")
SEEDS = (41242, 41243, 41244)
ROUTES = [ce.ROUTE_CLASSIC, ce.ROUTE_NON_CLASSIC]
MAIN = {
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1", "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O", "Caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "Lidocaine": "CCN(CC)CC(=O)Nc1c(C)cccc1C", "Benzocaine": "CCOC(=O)c1ccc(N)cc1",
    "Propranolol": "CC(C)NCC(O)COc1cccc2ccccc12", "Diclofenac": "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
}
EXTENDED = {
    "Naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O", "Ketoprofen": "CC(C(=O)O)c1cccc(c1)C(=O)c1ccccc1",
    "Indomethacin": "COc1ccc2n(C(=O)c3ccc(Cl)cc3)c(C)c(CC(=O)O)c2c1", "Salicylic acid": "OC(=O)c1ccccc1O",
    "Phenacetin": "CCOc1ccc(NC(C)=O)cc1", "Procaine": "CCN(CC)CCOC(=O)c1ccc(N)cc1",
    "Nicotinamide": "NC(=O)c1cccnc1", "Theophylline": "Cn1c(=O)c2[nH]cnc2n(C)c1=O",
    "Metoprolol": "COCCc1ccc(OCC(O)CNC(C)C)cc1", "Atenolol": "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",
    "Warfarin": "CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O", "Diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12",
    "Fluoxetine": "CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1", "Celecoxib": "Cc1ccc(cc1)-c1cc(nn1-c1ccc(cc1)S(N)(=O)=O)C(F)(F)F",
    "Haloperidol": "OC1(CCN(CCCC(=O)c2ccc(F)cc2)CC1)c1ccc(Cl)cc1", "Metoclopramide": "CCN(CC)CCNC(=O)c1cc(Cl)c(N)cc1OC",
    "Tramadol": "COc1cccc(c1)C1(O)CCCCC1CN(C)C", "Phenytoin": "O=C1NC(=O)C(N1)(c1ccccc1)c1ccccc1",
    "Carbamazepine": "NC(=O)N1c2ccccc2C=Cc2ccccc12", "Ciprofloxacin": "OC(=O)C1=CN(C2CC2)c2cc(N3CCNCC3)c(F)cc2C1=O",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate_panel(panel: dict[str, str], target: int):
    records, runs = [], []
    for name, smiles in panel.items():
        for seed in SEEDS:
            log = []
            products = ce.generate_analogs(
                smiles, target, use_veber=True, selected_routes=ROUTES,
                log_callback=log.append, random_seed=seed,
                max_sa_score=6.0, allow_multistep=False,
            )
            for item in products:
                records.append({**item, "Prototype": name, "Seed": seed})
            matches = [re.search(r"Concluído: (\d+)/(\d+) análogos em (\d+) tentativas", line) for line in log]
            match = next((item for item in matches if item), None)
            runs.append({"Prototype": name, "Seed": seed, "Target": target,
                         "Returned": len(products), "Attempts": int(match.group(3)) if match else None,
                         "DistinctWithinRun": len({x["SMILES"] for x in products}),
                         "ExplorationEnabled": False})
    return records, runs


def write_panel(panel: dict[str, str], target: int, prefix: str, out: Path):
    records, runs = generate_panel(panel, target)
    pd.DataFrame(records).to_csv(out / f"records_{prefix}.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(runs).to_csv(out / f"run_metadata_{prefix}.csv", index=False, encoding="utf-8-sig")
    return records, runs


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    main_records, main_runs = write_panel(MAIN, 30, "main_benchmark", out)
    ext_records, ext_runs = write_panel(EXTENDED, 30, "extended_panel", out)
    direct_rows = []
    for name, smiles in MAIN.items():
        products = ce.generate_analogs(
            smiles, 200, use_veber=True, selected_routes=ROUTES,
            random_seed=SEEDS[0], max_sa_score=6.0, allow_multistep=False,
        )
        direct_rows.append({"Prototype": name, "DirectEligible": len(products),
                            "Target": 200, "ExplorationEnabled": False})
    pd.DataFrame(direct_rows).to_csv(out / "direct_enumeration_yield.csv", index=False)
    metadata = {
        "software": f"Analog-San {ve.APP_VERSION} (18 registered operations)", "validation": ve.VALIDATION_VERSION,
        "evidence_catalog": "empty; no reaction records catalogued; no route verification claimed",
        "Python": sys.version.split()[0], "RDKit": rdBase.rdkitVersion,
        "pandas": pd.__version__, "seeds": list(SEEDS), "target_per_run": 30,
        "main_prototypes": len(MAIN), "extended_prototypes": len(EXTENDED),
        "filters": {"Lipinski": "zero violations", "Veber": True,
                    "SA_score_max": 6.0, "high_risk_reactive_alerts": "excluded"},
        "multistep_exploration": False,
        "seed_interpretation": "Direct enumeration is deterministic; repeated seeds are checks of invariance, not independent samples.",
        "sha256": {},
    }
    for path in sorted(out.glob("*.csv")):
        metadata["sha256"][path.name] = sha256(path)
    (out / "metadata_main_benchmark.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"main={len(main_records)} records/{len(main_runs)} runs; extended={len(ext_records)} records/{len(ext_runs)} runs")


if __name__ == "__main__":
    main()
