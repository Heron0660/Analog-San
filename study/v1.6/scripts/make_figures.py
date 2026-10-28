"""Create publication figures from the current Analog-San rebenchmark outputs."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ORDER = ["Paracetamol", "Aspirin", "Ibuprofen", "Caffeine", "Lidocaine", "Benzocaine", "Propranolol", "Diclofenac"]
SMILES = {
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1", "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O", "Caffeine": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
    "Lidocaine": "CCN(CC)CC(=O)Nc1c(C)cccc1C", "Benzocaine": "CCOC(=O)c1ccc(N)cc1",
    "Propranolol": "CC(C)NCC(O)COc1cccc2ccccc12", "Diclofenac": "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
}
COLORS = dict(zip(ORDER, ["#7FA6CC", "#E6C48C", "#93BFA9", "#D9A9B3", "#DDA68B", "#B9AED6", "#A3ADB8", "#C4C892"]))
INK = "#35404A"


def save(fig, name: str, out: Path) -> None:
    fig.savefig(out / f"Figure_{name}.png", dpi=400, facecolor="white", bbox_inches="tight")
    fig.savefig(out / f"Figure_{name}.svg", facecolor="white", bbox_inches="tight")
    plt.close(fig)


def workflow(out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 2.5))
    ax.set_xlim(0, 100); ax.set_ylim(0, 34); ax.axis("off")
    steps = [
        ("Input", "prototype\nSMILES"), ("Enumerate", "18 direct\ngraph edits"),
        ("Sanitize", "RDKit validity +\ncanonical deduplication"),
        ("Triage", "Lipinski + Veber\nSA ≤ 6; reactive-alert exclusion"),
        ("Prioritize", "manual edit tier +\nheuristic metrics"),
        ("Export", "CSV / SMI / MOL2\n+ reports"),
    ]
    width, gap = 14.8, 2.0
    for i, (title, body) in enumerate(steps):
        x = 0.5 + i * (width + gap)
        ax.add_patch(FancyBboxPatch((x, 10), width, 19, boxstyle="round,pad=0.25,rounding_size=1.0", fc="#F0F4F6", ec="#738A99", lw=1.0))
        ax.text(x + width / 2, 24, title, ha="center", va="center", fontsize=8, weight="bold", color=INK)
        ax.text(x + width / 2, 16, body, ha="center", va="center", fontsize=6.5, color=INK, linespacing=1.3)
        if i < len(steps) - 1:
            ax.add_patch(FancyArrowPatch((x + width + 0.2, 19.5), (x + width + gap - 0.2, 19.5), arrowstyle="-|>", mutation_scale=8, color="#738A99", lw=1.0))
    ax.text(50, 4.0, "Evaluated protocol: multi-step exploration disabled. SA score and transformation tiers are heuristics, not route predictions.", ha="center", fontsize=7, color="#5D6872")
    save(fig, "1", out)


def create(data_dir: Path, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    main = pd.read_csv(data_dir / "records_main_benchmark.csv")
    unique = main.drop_duplicates(["Prototype", "SMILES"])
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})
    workflow(out)

    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    values = [unique.loc[unique.Prototype == p, "ScoreGlobal"].dropna().to_numpy() for p in ORDER]
    boxes = ax.boxplot(values, patch_artist=True, widths=0.58, showfliers=False, medianprops={"color": INK, "linewidth": 1.4})
    for box, proto in zip(boxes["boxes"], ORDER):
        box.set_facecolor(COLORS[proto]); box.set_alpha(0.8)
    ax.set_xticks(range(1, len(ORDER) + 1)); ax.set_xticklabels(ORDER, rotation=24, ha="right")
    ax.set_ylabel("Implementation-defined ScoreGlobal (0–100)"); ax.set_xlabel("Prototype")
    ax.grid(axis="y", color="#E1E5E8", linewidth=0.6); ax.set_axisbelow(True)
    save(fig, "2", out)

    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    markers = dict(zip(ORDER, ["o", "s", "^", "D", "v", "P", "X", "h"]))
    for proto in ORDER:
        group = unique[unique.Prototype == proto]
        ax.scatter(group.MW, group.LogP, s=31, marker=markers[proto], color=COLORS[proto], edgecolor=INK, linewidth=0.35, label=proto, alpha=0.85)
        parent = Chem.MolFromSmiles(SMILES[proto])
        ax.scatter([Descriptors.MolWt(parent)], [Descriptors.MolLogP(parent)], marker="*", s=135, color=COLORS[proto], edgecolor=INK, linewidth=0.7, zorder=5)
    ax.set_xlabel("Molecular weight (Da)"); ax.set_ylabel("RDKit Crippen MolLogP")
    ax.grid(color="#E1E5E8", linewidth=0.6); ax.set_axisbelow(True)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False, fontsize=7)
    save(fig, "3", out)

    tiers = sorted(unique.PrioridadeSintese.dropna().unique())
    table = pd.crosstab(unique.Prototype, unique.PrioridadeSintese).reindex(index=ORDER, columns=tiers, fill_value=0)
    fractions = table.div(table.sum(axis=1).replace(0, 1), axis=0)
    fig, ax = plt.subplots(figsize=(7.4, 3.8)); bottom = [0.0] * len(ORDER)
    palette = ["#DCE6F0", "#BED2E4", "#9DB9D4", "#7C9FC2", "#5C82AB", "#3F6590", "#2E4D70"]
    for i, tier in enumerate(tiers):
        vals = fractions[tier].to_numpy()
        ax.bar(range(len(ORDER)), vals, bottom=bottom, width=0.72, color=palette[min(i, len(palette)-1)], edgecolor="white", linewidth=0.6, label=f"Tier {int(tier)}")
        bottom = [a + b for a, b in zip(bottom, vals)]
    ax.set_xticks(range(len(ORDER))); ax.set_xticklabels(ORDER, rotation=24, ha="right")
    ax.set_ylim(0, 1); ax.set_ylabel("Fraction of unique direct analogues"); ax.set_xlabel("Prototype")
    ax.legend(title="Manual transformation tier", ncol=4, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), fontsize=8)
    save(fig, "4", out)

    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.7))
    examples = [
        (SMILES["Paracetamol"], "Paracetamol prototype"),
        ("CC(=O)N(C)c1ccc(O)cc1", "N-methyl edit"),
        ("COc1ccc(NC(C)=O)cc1", "O-methyl edit"),
    ]
    for axis, (smiles, title) in zip(axes, examples):
        axis.imshow(Draw.MolToImage(Chem.MolFromSmiles(smiles), size=(650, 420)))
        axis.axis("off"); axis.set_title(title, fontsize=9, color=INK)
    fig.subplots_adjust(wspace=0.03)
    save(fig, "5", out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--figure-dir", type=Path, required=True)
    args = parser.parse_args()
    create(args.data_dir.resolve(), args.figure_dir.resolve())


if __name__ == "__main__":
    main()
