from __future__ import annotations

import argparse
import hashlib
import json
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger, rdBase
from rdkit.Chem import rdFingerprintGenerator

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "study" / "v1.6" / "scripts"))
import chem_engine as ce
import validation_engine as ve
from run_crem_comparison import MAIN, EXTENDED

RDLogger.DisableLog("rdApp.*")
PROTOTYPES = {**MAIN, **EXTENDED}
SEED = 41242
FINGERPRINTS = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def metrics(items: list[dict]) -> dict:
    mols = [Chem.MolFromSmiles(item["SMILES"]) for item in items]
    fps = [FINGERPRINTS.GetFingerprint(mol) for mol in mols]
    nearest = []
    for index, fingerprint in enumerate(fps):
        similarities = [
            DataStructs.TanimotoSimilarity(fingerprint, other)
            for other_index, other in enumerate(fps) if index != other_index
        ]
        nearest.append(max(similarities) if similarities else 0.0)
    return {
        "N": len(items),
        "MeanNearestNeighborTanimoto": float(np.mean(nearest)) if nearest else None,
        "MeanParentTanimoto": float(np.mean([item["SimilaridadeTanimoto"] for item in items])) if items else None,
        "MeanSAScore": float(np.mean([item["SAScore"] for item in items])) if items else None,
        "TierCounts": {str(tier): sum(item["PrioridadeRegra"] == tier for item in items) for tier in range(1, 8)},
    }


def bootstrap_delta(values: list[float], seed: int = SEED, replicates: int = 20000) -> list[float]:
    rng = random.Random(seed)
    n = len(values)
    means = sorted(sum(rng.choice(values) for _ in range(n)) / n for _ in range(replicates))
    return [means[int(0.025 * replicates)], means[int(0.975 * replicates)]]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    per_prototype = []
    for name, smiles in PROTOTYPES.items():
        for target in (5, 10):
            common = dict(
                base_smiles=smiles,
                target_count=target,
                use_veber=True,
                selected_routes=[ce.ROUTE_CLASSIC, ce.ROUTE_NON_CLASSIC],
                random_seed=SEED,
                max_sa_score=6.0,
                allow_multistep=False,
            )
            for strategy in ("tier_priority", "tier_diverse"):
                candidates = ce.generate_analogs(**common, selection_strategy=strategy)
                row = {"Prototype": name, "Target": target, "Strategy": strategy, **metrics(candidates)}
                per_prototype.append(row)

    per_path = output_dir / "portfolio_selection_by_prototype.csv"
    pd.DataFrame(per_prototype).to_csv(per_path, index=False, encoding="utf-8-sig")
    summary = {"software": ve.APP_VERSION, "validation": ve.VALIDATION_VERSION,
               "RDKit": rdBase.rdkitVersion, "Python": sys.version.split()[0],
               "prototypes": len(PROTOTYPES), "seed": SEED,
               "filters": {"Lipinski": True, "Veber": True, "SA_score_max": 6.0,
                           "selected_reactive_alerts_excluded": True},
               "metric": "mean over candidates of each candidate's maximum Morgan radius-2, 2048-bit Tanimoto to any other selected candidate; lower indicates lower nearest-neighbour similarity",
               "bootstrap": {"unit": "prototype; paired by prototype", "replicates": 20000, "seed": SEED},
               "outputs_sha256": {per_path.name: sha256(per_path)}, "results": {}}
    data = pd.DataFrame(per_prototype)
    for target in (5, 10):
        subset = data[data.Target == target].pivot(index="Prototype", columns="Strategy")
        target_summary = {}
        for metric in ("MeanNearestNeighborTanimoto", "MeanParentTanimoto", "MeanSAScore"):
            delta = (subset[metric]["tier_diverse"] - subset[metric]["tier_priority"]).astype(float).tolist()
            target_summary[metric] = {"mean_paired_delta": float(np.mean(delta)),
                                      "bootstrap_95pct_CI": bootstrap_delta(delta),
                                      "improved_count": sum(value < 0 for value in delta),
                                      "worsened_count": sum(value > 0 for value in delta)}
        tiers = subset["TierCounts"]
        target_summary["tier_distribution_changed"] = int(sum(a != b for a, b in zip(
            tiers["tier_priority"], tiers["tier_diverse"])))
        summary["results"][str(target)] = target_summary
    summary_path = output_dir / "portfolio_selection_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "study" / "v1.7" / "outputs")
    run(parser.parse_args().output_dir.resolve())


if __name__ == "__main__":
    main()
