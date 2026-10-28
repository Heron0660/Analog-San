"""Summarize the paired CReM comparison and create Figure 6."""
from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "study" / "dados"
FIGURES = ROOT / "study" / "figures"


def main() -> None:
    crem_pool = pd.read_csv(DATA / "records_crem_comparison.csv")
    crem = crem_pool[crem_pool["SelectedFor100"] == 1].copy()
    analogsan = pd.read_csv(DATA / "records_analogsan_comparison.csv")
    metadata = pd.read_csv(DATA / "run_metadata_comparison.csv")
    records = pd.concat([crem, analogsan], ignore_index=True)
    records["Returned"] = 1

    summaries = []
    for tool, group in records.groupby("Tool", sort=False):
        runs = metadata[metadata["Tool"] == tool]
        summaries.append({
            "Tool": tool,
            "Runs": len(runs),
            "Returned records": len(group),
            "Pooled unique SMILES": group["SMILES"].nunique(),
            "Runs reaching 100": int((runs["Returned"] == 100).sum()),
            "Mean records per run": float(runs["Returned"].mean()),
            "Mean parent Tanimoto": float(group["SimilarityTanimoto"].mean()),
            "SD parent Tanimoto": float(group["SimilarityTanimoto"].std(ddof=1)),
            "Mean SA score": float(group["SAScore"].mean()),
            "SD SA score": float(group["SAScore"].std(ddof=1)),
            "Mean alert count": float(group["AlertCount"].mean()),
            "Mean molecular weight": float(group["MW"].mean()),
            "Mean MolLogP": float(group["LogP"].mean()),
        })
    pd.DataFrame(summaries).to_csv(DATA / "comparison_summary.csv", index=False, float_format="%.4f")

    by_prototype = records.groupby(["Prototype", "Tool"], as_index=False).agg(
        records=("SMILES", "size"), pooled_unique=("SMILES", "nunique"),
        mean_similarity=("SimilarityTanimoto", "mean"), mean_SA=("SAScore", "mean"),
        mean_alerts=("AlertCount", "mean"), mean_MW=("MW", "mean"), mean_LogP=("LogP", "mean"),
    )
    run_summary = metadata.groupby(["Prototype", "Tool"], as_index=False).agg(
        mean_returned=("Returned", "mean"), min_returned=("Returned", "min"),
        runs_at_target=("Returned", lambda x: int((x == 100).sum())),
    )
    by_prototype.merge(run_summary, on=["Prototype", "Tool"], how="left").to_csv(
        DATA / "comparison_summary_by_prototype.csv", index=False, float_format="%.4f"
    )

    paired = []
    for (prototype, seed), frame in records.groupby(["Prototype", "Seed"]):
        a = set(frame.loc[frame["Tool"] == "Analog-San", "SMILES"])
        c = set(frame.loc[frame["Tool"] == "CReM", "SMILES"])
        union = a | c
        paired.append({
            "Prototype": prototype, "Seed": seed, "Analog-San count": len(a),
            "CReM count": len(c), "Exact overlap": len(a & c),
            "Jaccard": len(a & c) / len(union) if union else 0.0,
        })
    pd.DataFrame(paired).to_csv(DATA / "comparison_pairwise_overlap.csv", index=False, float_format="%.4f")

    colors = {"Analog-San": "#2E7D5B", "CReM": "#7556A6"}
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.9), constrained_layout=True)
    panels = [
        (axes[0], "Returned structures per run", [metadata.loc[metadata.Tool == name, "Returned"] for name in ("Analog-San", "CReM")], (0, 105), ["Analog-San", "CReM"]),
        (axes[1], "Parent Morgan Tanimoto", [records.loc[records.Tool == name, "SimilarityTanimoto"] for name in ("Analog-San", "CReM")], (0, 1), ["Analog-San", "CReM"]),
        (axes[2], "RDKit synthetic-accessibility score", [records.loc[records.Tool == name, "SAScore"] for name in ("Analog-San", "CReM")], None, ["Analog-San", "CReM"]),
    ]
    for ax, title, values, limits, tick_labels in panels:
        boxes = ax.boxplot(values, tick_labels=tick_labels, patch_artist=True, showfliers=False)
        for box, name in zip(boxes["boxes"], ("Analog-San", "CReM")):
            box.set_facecolor(colors[name])
            box.set_alpha(0.78)
        for median in boxes["medians"]:
            median.set_color("#222222")
            median.set_linewidth(1.4)
        for position, (series, name) in enumerate(zip(values, ("Analog-San", "CReM")), start=1):
            x = [position + ((i % 9) - 4) * 0.012 for i in range(len(series))]
            ax.scatter(x, series, s=5 if len(series) > 100 else 18, color="#FFFFFF", edgecolor=colors[name], linewidth=0.35, alpha=0.62, zorder=3, rasterized=True)
        ax.set_title(title, fontsize=10)
        ax.grid(axis="y", color="#D8D8D8", linewidth=0.6, alpha=0.8)
        ax.set_axisbelow(True)
        ax.tick_params(axis="x", labelrotation=15)
        if limits:
            ax.set_ylim(*limits)
    axes[0].axhline(100, color="#555555", linestyle="--", linewidth=1)
    fig.suptitle("Paired analogue-generation benchmark across 28 prototypes", fontsize=12)
    fig.savefig(FIGURES / "figure_6_comparison.png", dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(FIGURES / "figure_6_comparison.svg", facecolor="white", bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
