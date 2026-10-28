"""Summarise the paired v1.4 comparison at the prototype level and plot Figure 6."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLORS = {"Analog-San": "#237A57", "CReM": "#7354A3"}
TOOLS = ("Analog-San", "CReM")


def summarize(data_dir: Path, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    crem_pool = pd.read_csv(data_dir / "records_crem_comparison.csv")
    analog_pool = pd.read_csv(data_dir / "records_analogsan_comparison.csv")
    runs = pd.read_csv(data_dir / "run_metadata_comparison.csv")
    crem = crem_pool[crem_pool.SelectedForTarget == 1].copy()
    analog = analog_pool.copy()
    records = pd.concat([analog, crem], ignore_index=True)

    per_prototype = []
    for prototype in runs.Prototype.drop_duplicates():
        for tool in TOOLS:
            run_group = runs[(runs.Prototype == prototype) & (runs.Tool == tool)]
            molecule_group = records[(records.Prototype == prototype) & (records.Tool == tool)]
            per_prototype.append({
                "Prototype": prototype,
                "Tool": tool,
                "MeanReturned": run_group.Returned.mean(),
                "MinReturned": run_group.Returned.min(),
                "MaxReturned": run_group.Returned.max(),
                "RunsAtTarget": int((run_group.Returned == run_group.Target).sum()),
                "PooledUnique": molecule_group.SMILES.nunique(),
                "MeanParentTanimoto": molecule_group.SimilarityTanimoto.mean(),
                "MeanSAScore": molecule_group.SAScore.mean(),
                "MeanAlertCount": molecule_group.AlertCount.mean(),
                "Records": len(molecule_group),
            })
    by_prototype = pd.DataFrame(per_prototype)
    by_prototype.to_csv(data_dir / "comparison_summary_by_prototype.csv", index=False, float_format="%.4f")

    overall = []
    for tool in TOOLS:
        group = records[records.Tool == tool]
        run_group = runs[runs.Tool == tool]
        overall.append({
            "Tool": tool,
            "Prototypes": run_group.Prototype.nunique(),
            "Runs": len(run_group),
            "ReturnedRecords": int(run_group.Returned.sum()),
            "PooledUniqueSMILES": group.SMILES.nunique(),
            "RunsAtTarget": int((run_group.Returned == run_group.Target).sum()),
            "MeanReturnedPerRun": run_group.Returned.mean(),
            "MeanParentTanimoto": group.SimilarityTanimoto.mean(),
            "MeanSAScore": group.SAScore.mean(),
            "MeanAlertCount": group.AlertCount.mean(),
            "MeanPooledUniquePerPrototype": by_prototype[by_prototype.Tool == tool].PooledUnique.mean(),
        })
    pd.DataFrame(overall).to_csv(data_dir / "comparison_summary.csv", index=False, float_format="%.4f")

    # One point per prototype: seeds are repeated computational conditions, not independent samples.
    pivot = by_prototype.pivot(index="Prototype", columns="Tool")
    paired = []
    metrics = {
        "Returned": "MeanReturned",
        "PooledUnique": "PooledUnique",
        "ParentTanimoto": "MeanParentTanimoto",
        "SAScore": "MeanSAScore",
    }
    for prototype in pivot.index:
        row = {"Prototype": prototype}
        for label, column in metrics.items():
            row[f"AnalogSan_{label}"] = pivot.loc[prototype, (column, "Analog-San")]
            row[f"CReM_{label}"] = pivot.loc[prototype, (column, "CReM")]
            row[f"Difference_{label}"] = row[f"AnalogSan_{label}"] - row[f"CReM_{label}"]
        paired.append(row)
    pd.DataFrame(paired).to_csv(data_dir / "comparison_paired_by_prototype.csv", index=False, float_format="%.4f")

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.fonttype": "none"})
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 6.0), constrained_layout=True)
    panels = [
        ("MeanReturned", "Returned structures per run", "Structures / run"),
        ("PooledUnique", "Distinct structures across 3 runs", "Unique SMILES / prototype"),
        ("MeanParentTanimoto", "Similarity to prototype", "Morgan Tanimoto"),
        ("MeanSAScore", "RDKit SA Score (heuristic; lower is lower)", "SA score"),
    ]
    for ax, (metric, title, ylabel) in zip(axes.flat, panels):
        data = [by_prototype.loc[by_prototype.Tool == tool, metric].dropna().to_numpy() for tool in TOOLS]
        boxes = ax.boxplot(data, tick_labels=TOOLS, patch_artist=True, showfliers=False, widths=0.55)
        for box, tool in zip(boxes["boxes"], TOOLS):
            box.set_facecolor(COLORS[tool])
            box.set_alpha(0.72)
        for med in boxes["medians"]:
            med.set_color("#20252A")
            med.set_linewidth(1.4)
        for pos, (vals, tool) in enumerate(zip(data, TOOLS), start=1):
            jitter = [pos + ((i % 7) - 3) * 0.014 for i in range(len(vals))]
            ax.scatter(jitter, vals, s=15, color="white", edgecolor=COLORS[tool],
                       linewidth=0.55, alpha=0.8, zorder=3)
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.grid(axis="y", color="#D9DEE3", linewidth=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", labelrotation=12)
    fig.suptitle("Matched-input descriptive comparison (28 prototypes)", fontsize=12)
    fig.savefig(figure_dir / "Figure_6.png", dpi=400, facecolor="white", bbox_inches="tight")
    fig.savefig(figure_dir / "Figure_6.svg", facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    args = parser.parse_args()
    summarize(args.data_dir.resolve(), args.figure_dir.resolve())


if __name__ == "__main__":
    main()
