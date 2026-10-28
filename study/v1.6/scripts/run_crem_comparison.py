"""Paired comparison of Analog-San v1.3 and CReM 0.3.2 on 28 prototypes."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import random
import sys
from pathlib import Path

from rdkit import Chem, RDLogger, rdBase
from crem import __version__ as crem_version
from crem.crem import grow_mol, mutate_mol

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
import chem_engine as ce
import validation_engine as ve

RDLogger.DisableLog("rdApp.*")

SEEDS = (41242, 41243, 41244)
TARGET = 100
MAX_REPLACEMENTS_PER_MODE = 1000
MAIN = {
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O",
    "Caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "Lidocaine": "CCN(CC)CC(=O)Nc1c(C)cccc1C",
    "Benzocaine": "CCOC(=O)c1ccc(N)cc1",
    "Propranolol": "CC(C)NCC(O)COc1cccc2ccccc12",
    "Diclofenac": "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
}
EXTENDED = {
    "Naproxen": "COc1ccc2cc(ccc2c1)C(C)C(=O)O",
    "Ketoprofen": "CC(C(=O)O)c1cccc(c1)C(=O)c1ccccc1",
    "Indomethacin": "COc1ccc2n(C(=O)c3ccc(Cl)cc3)c(C)c(CC(=O)O)c2c1",
    "Salicylic acid": "OC(=O)c1ccccc1O",
    "Phenacetin": "CCOc1ccc(NC(C)=O)cc1",
    "Procaine": "CCN(CC)CCOC(=O)c1ccc(N)cc1",
    "Nicotinamide": "NC(=O)c1cccnc1",
    "Theophylline": "Cn1c(=O)c2[nH]cnc2n(C)c1=O",
    "Metoprolol": "COCCc1ccc(OCC(O)CNC(C)C)cc1",
    "Atenolol": "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1",
    "Warfarin": "CC(=O)CC(c1ccccc1)C1=C(O)c2ccccc2OC1=O",
    "Diazepam": "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc12",
    "Fluoxetine": "CNCCC(Oc1ccc(cc1)C(F)(F)F)c1ccccc1",
    "Celecoxib": "Cc1ccc(cc1)-c1cc(nn1-c1ccc(cc1)S(N)(=O)=O)C(F)(F)F",
    "Haloperidol": "OC1(CCN(CCCC(=O)c2ccc(F)cc2)CC1)c1ccc(Cl)cc1",
    "Metoclopramide": "CCN(CC)CCNC(=O)c1cc(Cl)c(N)cc1OC",
    "Tramadol": "COc1cccc(c1)C1(O)CCCCC1CN(C)C",
    "Phenytoin": "O=C1NC(=O)C(N1)(c1ccccc1)c1ccccc1",
    "Carbamazepine": "NC(=O)N1c2ccccc2C=Cc2ccccc12",
    "Ciprofloxacin": "OC(=O)C1=CN(C2CC2)c2cc(N3CCNCC3)c(F)cc2C1=O",
}
PROTOTYPES = {**MAIN, **EXTENDED}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def describe(name: str, seed: int, tool: str, mol: Chem.Mol) -> dict:
    props = ce.calculate_properties(mol)
    base = Chem.MolFromSmiles(PROTOTYPES[name])
    return {
        "Prototype": name,
        "Seed": seed,
        "Tool": tool,
        "SMILES": Chem.MolToSmiles(mol, isomericSmiles=True),
        **props,
        "SimilarityTanimoto": ve.tanimoto_similarity(base, mol),
        "SAScore": ve.synthetic_accessibility_score(mol),
        "AlertCount": len(ve.chemical_alerts(mol)),
    }


def passes_filters(mol: Chem.Mol) -> bool:
    props = ce.calculate_properties(mol)
    return ce.lipinski_filter(props) and ce.veber_filter(props)


def collect_crem(smiles: str, db_path: str, seed: int) -> tuple[list[str], dict]:
    parent = Chem.MolFromSmiles(smiles)
    if parent is None:
        raise ValueError(f"Invalid prototype SMILES: {smiles}")
    raw_by_mode = {}
    invalid_by_mode = {}
    structures = set()
    for mode, generator, kwargs in (
        ("mutate", mutate_mol, {"radius": 3, "min_size": 0, "max_size": 1, "max_inc": 3}),
        ("grow", grow_mol, {"radius": 3, "min_atoms": 1, "max_atoms": 1}),
    ):
        products = list(generator(
            parent,
            db_name=db_path,
            max_replacements=MAX_REPLACEMENTS_PER_MODE,
            seed=seed,
            **kwargs,
        ))
        raw_by_mode[mode] = len(products)
        invalid_by_mode[mode] = 0
        for product in products:
            mol = Chem.MolFromSmiles(product)
            if mol is None:
                invalid_by_mode[mode] += 1
                continue
            canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
            if canonical != Chem.MolToSmiles(parent, isomericSmiles=True):
                structures.add(canonical)
    raw_by_mode["invalid"] = sum(invalid_by_mode.values())
    raw_by_mode["invalid_mutate"] = invalid_by_mode["mutate"]
    raw_by_mode["invalid_grow"] = invalid_by_mode["grow"]
    return sorted(structures), raw_by_mode


def write_rows(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(db_path: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    crem_records: list[dict] = []
    analogsan_records: list[dict] = []
    metadata: list[dict] = []

    for name, smiles in PROTOTYPES.items():
        for seed in SEEDS:
            candidates, mode_counts = collect_crem(smiles, str(db_path), seed)
            accepted = [
                canonical for canonical in candidates
                if (mol := Chem.MolFromSmiles(canonical)) is not None and passes_filters(mol)
            ]
            selected = random.Random(seed).sample(accepted, min(TARGET, len(accepted)))
            crem_records.extend(
                {
                    **describe(name, seed, "CReM", Chem.MolFromSmiles(smi)),
                    "SelectedFor100": int(smi in selected),
                }
                for smi in accepted
            )

            generated = ce.generate_analogs(
                smiles,
                TARGET,
                use_veber=True,
                selected_routes=[ce.ROUTE_CLASSIC, ce.ROUTE_NON_CLASSIC],
                random_seed=seed,
            )
            analogsan_records.extend(
                describe(name, seed, "Analog-San", Chem.MolFromSmiles(item["SMILES"]))
                for item in generated
            )
            metadata.extend([
                {
                    "Prototype": name, "Seed": seed, "Tool": "CReM",
                    "RawMutate": mode_counts["mutate"], "RawGrow": mode_counts["grow"],
                    "InvalidProducts": mode_counts["invalid"],
                    "InvalidMutate": mode_counts["invalid_mutate"],
                    "InvalidGrow": mode_counts["invalid_grow"],
                    "UniqueCandidates": len(candidates), "PassBothFilters": len(accepted),
                    "Returned": len(selected), "Target": TARGET,
                },
                {
                    "Prototype": name, "Seed": seed, "Tool": "Analog-San",
                    "RawMutate": "NA", "RawGrow": "NA", "InvalidProducts": "NA",
                    "InvalidMutate": "NA", "InvalidGrow": "NA", "UniqueCandidates": "NA",
                    "PassBothFilters": "NA", "Returned": len(generated), "Target": TARGET,
                },
            ])

    write_rows(output_dir / "records_crem_comparison.csv", crem_records)
    write_rows(output_dir / "records_analogsan_comparison.csv", analogsan_records)
    write_rows(output_dir / "run_metadata_comparison.csv", metadata)
    metadata_doc = {
        "protocol": "Paired 28-prototype comparison; 3 fixed seeds; target 100; strict Lipinski and Veber filters.",
        "CReM": crem_version,
        "RDKit": rdBase.rdkitVersion,
        "Python": sys.version.split()[0],
        "database": "CReM ChEMBL33 SA<=2, alert-filtered, max fragment size 15; Zenodo 10.5281/zenodo.16909329; file chembl33_sa2_f5.db.gz",
        "database_sha256_gz": sha256(db_path.with_suffix(db_path.suffix + ".gz")),
        "parameters": {
            "mutate": {"radius": 3, "min_size": 0, "max_size": 1, "max_inc": 3, "max_replacements_per_seed": MAX_REPLACEMENTS_PER_MODE},
            "grow": {"radius": 3, "min_atoms": 1, "max_atoms": 1, "max_replacements_per_seed": MAX_REPLACEMENTS_PER_MODE},
            "CReM_selection": "Uniform seeded sample without replacement from unique structures passing both filters, capped at 100 per run.",
            "Analog-San": "Official v1.3 generator; both routes; Lipinski and Veber active; target 100.",
            "seeds": list(SEEDS),
        },
        "database_url": "https://zenodo.org/records/16909329/files/chembl33_sa2_f5.db.gz?download=1",
        "outputs_sha256": {},
    }
    for filename in ("records_crem_comparison.csv", "records_analogsan_comparison.csv", "run_metadata_comparison.csv"):
        metadata_doc["outputs_sha256"][filename] = sha256(output_dir / filename)
    (output_dir / "metadata_comparison.json").write_text(json.dumps(metadata_doc, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", type=Path, required=True, help="Path to the extracted CReM SQLite database")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "study" / "dados")
    args = parser.parse_args()
    if not args.db.is_file():
        parser.error(f"CReM database not found: {args.db}")
    run(args.db.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()
