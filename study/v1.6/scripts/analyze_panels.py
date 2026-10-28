"""Compute prototype-level descriptive tables for the current deterministic panels."""
from __future__ import annotations

import argparse
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")
FINGERPRINT = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)


def diversity(smiles: list[str]) -> float:
    fps = [FINGERPRINT.GetFingerprint(Chem.MolFromSmiles(smi)) for smi in smiles]
    pairs = [1.0 - DataStructs.TanimotoSimilarity(a, b) for a, b in itertools.combinations(fps, 2)]
    return float(np.mean(pairs)) if pairs else 0.0


def summarize_panel(records_path: Path, runs_path: Path, direct_path: Path | None = None):
    records = pd.read_csv(records_path)
    runs = pd.read_csv(runs_path)
    rows, diversity_rows = [], []
    for prototype, group in records.groupby("Prototype", sort=False):
        unique = group.drop_duplicates("SMILES")
        by_seed = [set(g.SMILES) for _, g in group.groupby("Seed")]
        union = set.union(*by_seed) if by_seed else set()
        shared = set.intersection(*by_seed) if by_seed else set()
        run_group = runs[runs.Prototype == prototype]
        rows.append({
            "Prototype": prototype, "Records": len(group), "DistinctAcrossSeeds": len(union),
            "SeedInvariant": len(union) == len(shared),
            "MeanScore": unique.ScoreGlobal.mean(), "SDScore": unique.ScoreGlobal.std(ddof=1),
            "MeanTanimoto": unique.SimilaridadeTanimoto.mean(), "MeanSA": unique.SAScore.mean(),
            "AlertFreePct": 100.0 * (unique.QtdAlertas == 0).mean(),
            "MeanReturnedPerRun": run_group.Returned.mean(),
            "RunsAtTarget": int((run_group.Returned == run_group.Target).sum()),
        })
        jaccards = [len(a & b) / len(a | b) if a | b else 1.0 for a, b in itertools.combinations(by_seed, 2)]
        scaffolds = {MurckoScaffold.MurckoScaffoldSmiles(mol=Chem.MolFromSmiles(smi)) for smi in union}
        diversity_rows.append({
            "Prototype": prototype,
            "DirectPerRun": float(group.groupby("Seed").size().mean()),
            "ExplorationPerRun": 0.0,
            "DistinctAcrossSeeds": len(union), "SharedAcrossSeeds": len(shared),
            "MeanPairwiseJaccard": float(np.mean(jaccards)) if jaccards else 1.0,
            "InternalDiversity": float(np.mean([diversity(sorted(values)) for values in by_seed])) if by_seed else 0.0,
            "MurckoScaffolds": len(scaffolds),
        })
    return pd.DataFrame(rows), pd.DataFrame(diversity_rows), records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    args = parser.parse_args()
    root = args.data_dir.resolve()
    for panel in ("main_benchmark", "extended_panel"):
        summary, diversity_table, records = summarize_panel(
            root / f"records_{panel}.csv", root / f"run_metadata_{panel}.csv"
        )
        summary.to_csv(root / f"summary_{panel}.csv", index=False, float_format="%.4f")
        diversity_table.to_csv(root / f"diversity_{panel}.csv", index=False, float_format="%.4f")
        print(panel, "records", len(records), "unique", records.SMILES.nunique(),
              "mean score", records.ScoreGlobal.mean(), "SD", records.ScoreGlobal.std(),
              "mean similarity", records.SimilaridadeTanimoto.mean(),
              "mean SA", records.SAScore.mean(), "alert-free", int((records.QtdAlertas == 0).sum()))
        print(summary.round(3).to_string(index=False))
    print("direct yields")
    print(pd.read_csv(root / "direct_enumeration_yield.csv").to_string(index=False))


if __name__ == "__main__":
    main()
