"""
=============================================================================
  VALIDAÇÃO CIENTÍFICA - Analog-San
  Módulo: validation_engine.py
=============================================================================
  Calcula métricas heurísticas de priorização, SA Score estrutural,
  similaridade ao composto-protótipo, alertas químicos e gera relatórios com
  figuras simples para uso em artigo científico.
=============================================================================
"""

import csv
import hashlib
import json
import math
import os
from datetime import datetime
from pathlib import Path
from statistics import mean
from typing import Dict, List, Optional, Tuple

from PIL import Image, ImageDraw, ImageFont
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, FilterCatalog, rdFingerprintGenerator, rdMolDescriptors

try:
    from rdkit.Contrib.SA_Score import sascorer
except Exception:
    sascorer = None


VALIDATION_VERSION = "1.4"
APP_VERSION = "1.7.1"
DEFAULT_RANDOM_SEED = 0xA11A

RULE_PRIORITY_SCORE = {
    "muito alta": 100.0,
    "alta": 82.0,
    "média": 64.0,
    "media": 64.0,
    "baixa": 42.0,
}

SMARTS_ALERTS = [
    ("Nitro potencialmente crítico", "[N+](=O)[O-]"),
    ("Cloreto de acila reativo", "C(=O)Cl"),
    ("Cloreto de sulfonila reativo", "S(=O)(=O)Cl"),
    ("Isocianato reativo", "N=C=O"),
    ("Aldeído reativo", "[CX3H1](=O)[#6]"),
    ("Epóxido tensionado", "C1OC1"),
    ("Azida/diazo energético", "[N-]=[N+]=N"),
    ("Michael acceptor", "[C,c]=[C,c]-[C;X3](=O)"),
]

_SMARTS_ALERT_MOLS = [
    (name, Chem.MolFromSmarts(smarts)) for name, smarts in SMARTS_ALERTS
]
_MORGAN_GENERATOR = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
_FILTER_CATALOG = None


# ---------------------------------------------------------------------------
#  MÉTRICAS MOLECULARES
# ---------------------------------------------------------------------------

def _get_filter_catalog():
    """Inicializa PAINS/Brenk sob demanda para reduzir custo de abertura."""
    global _FILTER_CATALOG
    if _FILTER_CATALOG is not None:
        return _FILTER_CATALOG

    params = FilterCatalog.FilterCatalogParams()
    catalogs = FilterCatalog.FilterCatalogParams.FilterCatalogs
    for catalog_name in ("PAINS_A", "PAINS_B", "PAINS_C", "BRENK"):
        if hasattr(catalogs, catalog_name):
            params.AddCatalog(getattr(catalogs, catalog_name))
    _FILTER_CATALOG = FilterCatalog.FilterCatalog(params)
    return _FILTER_CATALOG


def _clip(value: float, low: float = 0.0, high: float = 100.0) -> float:
    return max(low, min(high, value))


def _range_score(value: float, low: float, high: float, tolerance: float) -> float:
    """Score 0-100: máximo dentro da faixa, queda linear fora dela."""
    if low <= value <= high:
        return 100.0
    if value < low:
        return _clip(100.0 - ((low - value) / tolerance) * 100.0)
    return _clip(100.0 - ((value - high) / tolerance) * 100.0)


def _lipinski_violations(props: Dict) -> List[str]:
    violations = []
    if props.get("MW", 0) > 500:
        violations.append("MW>500")
    if props.get("LogP", 0) > 5:
        violations.append("LogP>5")
    if props.get("HBD", 0) > 5:
        violations.append("HBD>5")
    if props.get("HBA", 0) > 10:
        violations.append("HBA>10")
    return violations


def _veber_violations(props: Dict) -> List[str]:
    violations = []
    if props.get("RotBonds", 0) > 10:
        violations.append("RotBonds>10")
    if props.get("TPSA", 0) > 140:
        violations.append("TPSA>140")
    return violations


def _fingerprint(mol: Chem.Mol):
    return _MORGAN_GENERATOR.GetFingerprint(mol)


def tanimoto_similarity(base_mol: Chem.Mol, analog_mol: Chem.Mol) -> float:
    try:
        return round(float(DataStructs.TanimotoSimilarity(
            _fingerprint(base_mol), _fingerprint(analog_mol)
        )), 4)
    except Exception:
        return 0.0


def synthetic_accessibility_score(mol: Chem.Mol) -> float:
    """
    Retorna SA Score estimado. Na escala clássica, 1 é mais fácil e 10 mais
    difícil. Usa o sascorer do RDKit quando disponível e fallback heurístico.
    """
    if sascorer is not None:
        try:
            return round(float(sascorer.calculateScore(mol)), 2)
        except Exception:
            pass

    hac = max(1, mol.GetNumHeavyAtoms())
    ring_count = rdMolDescriptors.CalcNumRings(mol)
    chiral = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
    rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    hetero = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in (1, 6))
    raw = 1.0 + hac / 18.0 + ring_count * 0.35 + chiral * 0.45 + rot * 0.08 + hetero * 0.04
    return round(_clip(raw, 1.0, 10.0), 2)


def chemical_alerts(mol: Chem.Mol) -> List[str]:
    alerts: List[str] = []

    try:
        matches = _get_filter_catalog().GetMatches(mol)
        for entry in matches:
            desc = entry.GetDescription()
            if desc and desc not in alerts:
                alerts.append(desc)
    except Exception:
        pass

    for name, patt in _SMARTS_ALERT_MOLS:
        if patt is not None and mol.HasSubstructMatch(patt) and name not in alerts:
            alerts.append(name)

    halogens = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() in (9, 17, 35, 53))
    if halogens > 5:
        alerts.append("Excesso de halogênios")

    return alerts


def score_candidate(
    props: Dict,
    similarity: float,
    rule_priority: int,
    rule_tier: str,
    sa_score: float,
    alerts: List[str],
    use_veber: bool = True,
) -> Dict:
    lip_v = _lipinski_violations(props)
    veb_v = _veber_violations(props) if use_veber else []

    lipinski_score = _clip(100.0 - len(lip_v) * 25.0)
    veber_score = _clip(100.0 - len(veb_v) * 35.0)
    druglike_score = (lipinski_score * 0.65) + (veber_score * 0.35)

    physicochemical_score = mean([
        _range_score(float(props.get("MW", 0)), 180, 500, 180),
        _range_score(float(props.get("LogP", 0)), 0.0, 4.5, 2.0),
        _range_score(float(props.get("TPSA", 0)), 20, 120, 80),
    ])
    similarity_score = _clip(similarity * 100.0)
    rule_priority_score = RULE_PRIORITY_SCORE.get(
        str(rule_tier).lower(), _clip(108.0 - rule_priority * 11.0)
    )
    sa_component = _clip(100.0 - ((sa_score - 1.0) / 9.0) * 100.0)
    alert_score = _clip(100.0 - len(alerts) * 30.0)

    global_score = (
        druglike_score * 0.18
        + physicochemical_score * 0.18
        + similarity_score * 0.18
        + rule_priority_score * 0.18
        + sa_component * 0.14
        + alert_score * 0.14
    )

    if global_score >= 80:
        category = "prioritário"
    elif global_score >= 65:
        category = "promissor"
    elif global_score >= 50:
        category = "moderado"
    else:
        category = "baixo"

    return {
        "ScoreGlobal": round(global_score, 2),
        "CategoriaScore": category,
        "LipinskiOK": "sim" if not lip_v else "não",
        "VeberOK": "sim" if not veb_v else ("n/a" if not use_veber else "não"),
        "RegrasVioladas": "; ".join(lip_v + veb_v) if lip_v or veb_v else "nenhuma",
    }


def evaluate_candidate(
    base_mol: Chem.Mol,
    analog_mol: Chem.Mol,
    props: Dict,
    rule_priority: int,
    rule_tier: str,
    use_veber: bool = True,
) -> Dict:
    """Calcula validação completa para um único análogo."""
    similarity = tanimoto_similarity(base_mol, analog_mol)
    sa_score = synthetic_accessibility_score(analog_mol)
    alerts = chemical_alerts(analog_mol)
    score = score_candidate(
        props=props,
        similarity=similarity,
        rule_priority=rule_priority,
        rule_tier=rule_tier,
        sa_score=sa_score,
        alerts=alerts,
        use_veber=use_veber,
    )
    score.update({
        "SimilaridadeTanimoto": similarity,
        "SAScore": sa_score,
        "QtdAlertas": len(alerts),
        "AlertasQuimicos": "; ".join(alerts) if alerts else "nenhum",
    })
    return score


def sort_by_scientific_score(analogs: List[Dict]) -> List[Dict]:
    """Ordena por frentes de Pareto quando presentes; caso contrário, pelo score legado."""
    if any("FrentePareto" in analog for analog in analogs):
        return sorted(
            analogs,
            key=lambda a: (
                int(a.get("FrentePareto", 10**9)),
                float(a.get("SAScore", 10)),
                -float(a.get("SimilaridadeTanimoto", 0)),
                int(a.get("QtdAlertas", 10**6)),
                a.get("SMILES", ""),
            ),
        )
    return sorted(
        analogs,
        key=lambda a: (
            -float(a.get("ScoreGlobal", 0)),
            int(a.get("PrioridadeRegra", a.get("PrioridadeSintese", 99))),
            float(a.get("SAScore", 10)),
            -float(a.get("SimilaridadeTanimoto", 0)),
            a.get("SMILES", ""),
        ),
    )


OBJECTIVE_PROPERTIES = {
    "LogP": "LogP", "MW": "MW", "TPSA": "TPSA",
    "HBD": "HBD", "HBA": "HBA", "RotBonds": "RotBonds",
}
OBJECTIVE_DIRECTIONS = {"increase", "decrease", "maintain"}


def rank_by_objectives(
    analogs: List[Dict], base_properties: Dict, objective_profile: List[Dict[str, str]]
) -> List[Dict]:
    """Anota frentes não dominadas usando preferências do usuário e métricas auxiliares."""
    if not objective_profile:
        return analogs
    if len(objective_profile) > 2:
        raise ValueError("Defina no máximo dois objetivos moleculares.")

    objectives = []
    seen = set()
    for item in objective_profile:
        prop = item.get("property")
        direction = item.get("direction")
        if prop not in OBJECTIVE_PROPERTIES or direction not in OBJECTIVE_DIRECTIONS:
            raise ValueError(f"Objetivo inválido: {item!r}")
        if prop in seen:
            raise ValueError(f"A propriedade {prop} foi selecionada mais de uma vez.")
        seen.add(prop)
        objectives.append((prop, direction))

    def utility_vector(analog):
        values = []
        deltas = {}
        for prop, direction in objectives:
            delta = float(analog[prop]) - float(base_properties[prop])
            deltas[prop] = delta
            values.append(delta if direction == "increase" else -delta if direction == "decrease" else -abs(delta))
        values.extend((
            float(analog.get("SimilaridadeTanimoto", 0)),
            -float(analog.get("SAScore", 10)),
            -float(analog.get("QtdAlertas", 0)),
        ))
        analog["DeltaObjetivos"] = deltas
        analog["Objetivos"] = "; ".join(
            f"{prop} {direction}: {deltas[prop]:+.2f}" for prop, direction in objectives
        )
        return values

    vectors = [utility_vector(a) for a in analogs]
    dominates = [[] for _ in analogs]
    dominated_count = [0] * len(analogs)
    fronts = [[]]
    for i, left in enumerate(vectors):
        for j, right in enumerate(vectors):
            if i == j:
                continue
            if all(x >= y for x, y in zip(left, right)) and any(x > y for x, y in zip(left, right)):
                dominates[i].append(j)
            elif all(y >= x for x, y in zip(left, right)) and any(y > x for x, y in zip(left, right)):
                dominated_count[i] += 1
        if dominated_count[i] == 0:
            analogs[i]["FrentePareto"] = 1
            fronts[0].append(i)

    rank = 1
    while fronts[-1]:
        next_front = []
        for i in fronts[-1]:
            for j in dominates[i]:
                dominated_count[j] -= 1
                if dominated_count[j] == 0:
                    analogs[j]["FrentePareto"] = rank + 1
                    next_front.append(j)
        rank += 1
        fronts.append(next_front)
    return sort_by_scientific_score(analogs)


# ---------------------------------------------------------------------------
#  FIGURAS E RELATÓRIO
# ---------------------------------------------------------------------------

def _font(size: int, bold: bool = False):
    candidates = [
        "arialbd.ttf" if bold else "arial.ttf",
        "segoeuib.ttf" if bold else "segoeui.ttf",
    ]
    for name in candidates:
        try:
            return ImageFont.truetype(name, size)
        except Exception:
            continue
    return ImageFont.load_default()


def _text(draw: ImageDraw.ImageDraw, xy, text: str, fill, font):
    draw.text(xy, text, fill=fill, font=font)


def _new_canvas(width: int = 1400, height: int = 900) -> Tuple[Image.Image, ImageDraw.ImageDraw]:
    img = Image.new("RGB", (width, height), "#ffffff")
    draw = ImageDraw.Draw(img)
    return img, draw


def _save_png(img: Image.Image, path: Path) -> str:
    img.save(path, "PNG", dpi=(300, 300))
    return str(path)


def _svg_header(width: int, height: int, title: str) -> List[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        f'<text x="32" y="42" font-family="Arial" font-size="24" font-weight="700" fill="#111827">{_xml(title)}</text>',
    ]


def _xml(text: object) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def _write_svg(lines: List[str], path: Path) -> str:
    lines.append("</svg>")
    path.write_text("\n".join(lines), encoding="utf-8")
    return str(path)


def _axis(draw, box, x_label: str, y_label: str):
    x0, y0, x1, y1 = box
    draw.line((x0, y1, x1, y1), fill="#111827", width=2)
    draw.line((x0, y0, x0, y1), fill="#111827", width=2)
    _text(draw, (x0 + (x1 - x0) // 2 - 35, y1 + 28), x_label, "#374151", _font(18))
    _text(draw, (x0 - 58, y0 - 8), y_label, "#374151", _font(18))


def _histogram(draw, values: List[float], box, title: str, color: str):
    x0, y0, x1, y1 = box
    _text(draw, (x0, y0 - 38), title, "#111827", _font(20, True))
    _axis(draw, box, "valor", "n")
    if not values:
        return
    bins = min(10, max(4, int(math.sqrt(len(values)))))
    lo, hi = min(values), max(values)
    if math.isclose(lo, hi):
        lo -= 0.5
        hi += 0.5
    counts = [0] * bins
    for v in values:
        idx = min(bins - 1, int(((v - lo) / (hi - lo)) * bins))
        counts[idx] += 1
    max_count = max(counts) or 1
    bar_w = (x1 - x0 - 12) / bins
    for i, count in enumerate(counts):
        bh = (count / max_count) * (y1 - y0 - 20)
        bx0 = x0 + 7 + i * bar_w
        bx1 = bx0 + bar_w * 0.72
        by0 = y1 - bh
        draw.rectangle((bx0, by0, bx1, y1), fill=color)
    _text(draw, (x0, y1 + 4), f"{lo:.1f}", "#6b7280", _font(14))
    _text(draw, (x1 - 55, y1 + 4), f"{hi:.1f}", "#6b7280", _font(14))


def _scatter(draw, x_values, y_values, box, x_label: str, y_label: str, color: str):
    x0, y0, x1, y1 = box
    _axis(draw, box, x_label, y_label)
    if not x_values or not y_values:
        return
    xmin, xmax = min(x_values), max(x_values)
    ymin, ymax = min(y_values), max(y_values)
    if math.isclose(xmin, xmax):
        xmin -= 1
        xmax += 1
    if math.isclose(ymin, ymax):
        ymin -= 1
        ymax += 1
    for x, y in zip(x_values, y_values):
        px = x0 + ((x - xmin) / (xmax - xmin)) * (x1 - x0)
        py = y1 - ((y - ymin) / (ymax - ymin)) * (y1 - y0)
        draw.ellipse((px - 5, py - 5, px + 5, py + 5), fill=color, outline="#ffffff")
    _text(draw, (x0, y1 + 4), f"{xmin:.1f}", "#6b7280", _font(14))
    _text(draw, (x1 - 55, y1 + 4), f"{xmax:.1f}", "#6b7280", _font(14))
    _text(draw, (x0 - 48, y1 - 15), f"{ymin:.1f}", "#6b7280", _font(14))
    _text(draw, (x0 - 48, y0 - 5), f"{ymax:.1f}", "#6b7280", _font(14))


def _make_distribution_png(analogs: List[Dict], path: Path):
    img, draw = _new_canvas()
    _text(draw, (32, 28), "Distribuição de descritores físico-químicos", "#111827", _font(28, True))
    _histogram(draw, [float(a["MW"]) for a in analogs], (70, 150, 420, 720), "Peso molecular (MW)", "#2563eb")
    _histogram(draw, [float(a["LogP"]) for a in analogs], (525, 150, 875, 720), "LogP", "#059669")
    _histogram(draw, [float(a["TPSA"]) for a in analogs], (980, 150, 1330, 720), "TPSA", "#d97706")
    _text(draw, (70, 820), "Figura 1. Distribuições utilizadas para avaliar drug-likeness e espaço químico dos análogos.", "#374151", _font(18))
    return _save_png(img, path)


def _make_scatter_png(analogs: List[Dict], path: Path, x_key: str, y_key: str, title: str, x_label: str, y_label: str):
    img, draw = _new_canvas()
    _text(draw, (32, 28), title, "#111827", _font(28, True))
    _scatter(
        draw,
        [float(a.get(x_key, 0)) for a in analogs],
        [float(a.get(y_key, 0)) for a in analogs],
        (120, 120, 1260, 760),
        x_label,
        y_label,
        "#2563eb",
    )
    return _save_png(img, path)


def _make_ranking_png(analogs: List[Dict], path: Path):
    top = sort_by_scientific_score(analogs)[:10]
    pareto_mode = any("FrentePareto" in a for a in top)
    img, draw = _new_canvas()
    title = "10 candidatos em ordem Pareto" if pareto_mode else "10 candidatos por ScoreGlobal"
    _text(draw, (32, 28), title, "#111827", _font(28, True))
    x0, y0, x1, y1 = 360, 130, 1260, 760
    max_score = max([float(a.get("ScoreGlobal", 0)) for a in top] or [1])
    bar_h = 46
    for i, analog in enumerate(top):
        y = y0 + i * 62
        score = float(analog.get("ScoreGlobal", 0))
        width = (score / max_score) * (x1 - x0)
        draw.rectangle((x0, y, x0 + width, y + bar_h), fill="#1f6feb")
        _text(draw, (70, y + 11), analog.get("ID_Molecula", f"#{i+1}"), "#111827", _font(18, True))
        _text(draw, (170, y + 11), f"{score:.1f}", "#111827", _font(18))
        _text(draw, (x0 + 8, y + 11), analog.get("CategoriaScore", ""), "#ffffff", _font(18, True))
    caption = ("Ordem por frentes de Pareto; barras mostram ScoreGlobal auxiliar."
               if pareto_mode else "Ordenação por ScoreGlobal heurístico.")
    _text(draw, (70, 820), f"Figura 5. {caption}", "#374151", _font(18))
    return _save_png(img, path)


def _make_radar_png(analogs: List[Dict], path: Path):
    top = sort_by_scientific_score(analogs)[:5]
    pareto_mode = any("FrentePareto" in a for a in top)
    img, draw = _new_canvas()
    _text(draw, (32, 28), "Radar dos melhores análogos", "#111827", _font(28, True))
    center = (700, 455)
    radius = 260
    axes = ["Score", "Similaridade", "Prioridade regra", "SA", "Sem alertas"]
    angles = [(-90 + i * 360 / len(axes)) * math.pi / 180 for i in range(len(axes))]
    for r in (0.25, 0.5, 0.75, 1.0):
        pts = [(center[0] + math.cos(a) * radius * r, center[1] + math.sin(a) * radius * r) for a in angles]
        draw.line(pts + [pts[0]], fill="#d1d5db", width=1)
    for axis, angle in zip(axes, angles):
        end = (center[0] + math.cos(angle) * radius, center[1] + math.sin(angle) * radius)
        draw.line((center[0], center[1], end[0], end[1]), fill="#9ca3af", width=1)
        label = (center[0] + math.cos(angle) * (radius + 42), center[1] + math.sin(angle) * (radius + 42))
        _text(draw, (label[0] - 55, label[1] - 10), axis, "#374151", _font(16))
    colors = ["#2563eb", "#059669", "#d97706", "#dc2626", "#7c3aed"]
    for idx, analog in enumerate(top):
        vals = [
            float(analog.get("ScoreGlobal", 0)) / 100,
            float(analog.get("SimilaridadeTanimoto", 0)),
            RULE_PRIORITY_SCORE.get(str(analog.get("NivelRegraManual", "")).lower(), 60) / 100,
            _clip(100 - ((float(analog.get("SAScore", 10)) - 1) / 9) * 100) / 100,
            _clip(100 - float(analog.get("QtdAlertas", 0)) * 30) / 100,
        ]
        pts = [(center[0] + math.cos(a) * radius * v, center[1] + math.sin(a) * radius * v) for a, v in zip(angles, vals)]
        draw.line(pts + [pts[0]], fill=colors[idx % len(colors)], width=3)
        _text(draw, (100, 150 + idx * 32), analog.get("ID_Molecula", ""), colors[idx % len(colors)], _font(18, True))
    caption = ("Candidatos em ordem Pareto; Score é métrica auxiliar."
               if pareto_mode else "Candidatos ordenados pelo ScoreGlobal heurístico.")
    _text(draw, (80, 820), f"Figura 4. {caption}", "#374151", _font(18))
    return _save_png(img, path)


def _make_basic_svg(analogs: List[Dict], path: Path, title: str, x_key: str, y_key: str, x_label: str, y_label: str):
    width, height = 1000, 680
    lines = _svg_header(width, height, title)
    x0, y0, x1, y1 = 90, 90, 930, 590
    lines.append(f'<line x1="{x0}" y1="{y1}" x2="{x1}" y2="{y1}" stroke="#111827" stroke-width="2"/>')
    lines.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y1}" stroke="#111827" stroke-width="2"/>')
    xs = [float(a.get(x_key, 0)) for a in analogs]
    ys = [float(a.get(y_key, 0)) for a in analogs]
    if xs and ys:
        xmin, xmax = min(xs), max(xs)
        ymin, ymax = min(ys), max(ys)
        if math.isclose(xmin, xmax):
            xmin -= 1
            xmax += 1
        if math.isclose(ymin, ymax):
            ymin -= 1
            ymax += 1
        for x, y in zip(xs, ys):
            px = x0 + ((x - xmin) / (xmax - xmin)) * (x1 - x0)
            py = y1 - ((y - ymin) / (ymax - ymin)) * (y1 - y0)
            lines.append(f'<circle cx="{px:.1f}" cy="{py:.1f}" r="5" fill="#2563eb"/>')
    lines.append(f'<text x="{(x0+x1)/2-40:.1f}" y="{height-35}" font-family="Arial" font-size="18" fill="#374151">{_xml(x_label)}</text>')
    lines.append(f'<text x="18" y="{(y0+y1)/2:.1f}" font-family="Arial" font-size="18" fill="#374151">{_xml(y_label)}</text>')
    return _write_svg(lines, path)


def _make_distribution_svg(analogs: List[Dict], path: Path):
    width, height = 1000, 680
    lines = _svg_header(width, height, "Distribuição de descritores")
    keys = [("MW", "#2563eb", "MW"), ("LogP", "#059669", "LogP"), ("TPSA", "#d97706", "TPSA")]
    for panel, (key, color, label) in enumerate(keys):
        values = [float(a.get(key, 0)) for a in analogs]
        x0 = 80 + panel * 300
        y0, x1, y1 = 130, x0 + 235, 560
        lines.append(f'<text x="{x0}" y="112" font-family="Arial" font-size="18" font-weight="700">{label}</text>')
        lines.append(f'<line x1="{x0}" y1="{y1}" x2="{x1}" y2="{y1}" stroke="#111827" stroke-width="2"/>')
        lines.append(f'<line x1="{x0}" y1="{y0}" x2="{x0}" y2="{y1}" stroke="#111827" stroke-width="2"/>')
        if values:
            bins = min(10, max(4, int(math.sqrt(len(values)))))
            lo, hi = min(values), max(values)
            if math.isclose(lo, hi):
                lo -= 0.5
                hi += 0.5
            counts = [0] * bins
            for v in values:
                idx = min(bins - 1, int(((v - lo) / (hi - lo)) * bins))
                counts[idx] += 1
            max_count = max(counts) or 1
            bar_w = (x1 - x0 - 8) / bins
            for i, count in enumerate(counts):
                bh = (count / max_count) * (y1 - y0 - 18)
                bx = x0 + 6 + i * bar_w
                lines.append(f'<rect x="{bx:.1f}" y="{y1-bh:.1f}" width="{bar_w*0.72:.1f}" height="{bh:.1f}" fill="{color}"/>')
    return _write_svg(lines, path)


def _make_ranking_svg(analogs: List[Dict], path: Path):
    top = sort_by_scientific_score(analogs)[:10]
    width, height = 1000, 680
    title = "Ordem Pareto (barras: ScoreGlobal auxiliar)" if any("FrentePareto" in a for a in top) else "Ranking por ScoreGlobal heurístico"
    lines = _svg_header(width, height, title)
    max_score = max([float(a.get("ScoreGlobal", 0)) for a in top] or [1])
    for i, analog in enumerate(top):
        y = 95 + i * 53
        score = float(analog.get("ScoreGlobal", 0))
        w = (score / max_score) * 650
        lines.append(f'<text x="58" y="{y+24}" font-family="Arial" font-size="16" font-weight="700" fill="#111827">{_xml(analog.get("ID_Molecula", ""))}</text>')
        lines.append(f'<text x="172" y="{y+24}" font-family="Arial" font-size="16" fill="#111827">{score:.1f}</text>')
        lines.append(f'<rect x="250" y="{y}" width="{w:.1f}" height="34" fill="#1f6feb"/>')
        lines.append(f'<text x="260" y="{y+23}" font-family="Arial" font-size="14" font-weight="700" fill="#ffffff">{_xml(analog.get("CategoriaScore", ""))}</text>')
    return _write_svg(lines, path)


def _make_radar_svg(analogs: List[Dict], path: Path):
    top = sort_by_scientific_score(analogs)[:5]
    width, height = 1000, 680
    lines = _svg_header(width, height, "Radar dos melhores análogos")
    cx, cy, radius = 500, 365, 220
    axes = ["Score", "Similaridade", "Prioridade regra", "SA", "Sem alertas"]
    angles = [(-90 + i * 360 / len(axes)) * math.pi / 180 for i in range(len(axes))]
    for r in (0.25, 0.5, 0.75, 1.0):
        pts = " ".join(f"{cx + math.cos(a)*radius*r:.1f},{cy + math.sin(a)*radius*r:.1f}" for a in angles)
        lines.append(f'<polygon points="{pts}" fill="none" stroke="#d1d5db" stroke-width="1"/>')
    colors = ["#2563eb", "#059669", "#d97706", "#dc2626", "#7c3aed"]
    for idx, analog in enumerate(top):
        vals = [
            float(analog.get("ScoreGlobal", 0)) / 100,
            float(analog.get("SimilaridadeTanimoto", 0)),
            RULE_PRIORITY_SCORE.get(str(analog.get("NivelRegraManual", "")).lower(), 60) / 100,
            _clip(100 - ((float(analog.get("SAScore", 10)) - 1) / 9) * 100) / 100,
            _clip(100 - float(analog.get("QtdAlertas", 0)) * 30) / 100,
        ]
        pts = " ".join(f"{cx + math.cos(a)*radius*v:.1f},{cy + math.sin(a)*radius*v:.1f}" for a, v in zip(angles, vals))
        lines.append(f'<polygon points="{pts}" fill="none" stroke="{colors[idx % len(colors)]}" stroke-width="3"/>')
        lines.append(f'<text x="70" y="{120 + idx*28}" font-family="Arial" font-size="16" font-weight="700" fill="{colors[idx % len(colors)]}">{_xml(analog.get("ID_Molecula", ""))}</text>')
    for axis, angle in zip(axes, angles):
        lx = cx + math.cos(angle) * (radius + 38)
        ly = cy + math.sin(angle) * (radius + 38)
        lines.append(f'<text x="{lx-42:.1f}" y="{ly:.1f}" font-family="Arial" font-size="15" fill="#374151">{_xml(axis)}</text>')
    return _write_svg(lines, path)


def _result_hash(analogs: List[Dict]) -> str:
    payload = json.dumps(analogs, sort_keys=True, ensure_ascii=False).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def export_validated_csv(analogs: List[Dict], filepath: Path) -> str:
    if not analogs:
        raise ValueError("Nenhum análogo para exportar.")
    fieldnames = sorted({key for analog in analogs for key in analog.keys()})
    with filepath.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(analogs)
    return str(filepath)


def generate_scientific_report(
    analogs: List[Dict],
    base_smiles: str,
    output_dir: str,
    params: Optional[Dict] = None,
) -> Dict:
    """Gera relatório, CSV validado, metadados e figuras PNG/SVG."""
    if not analogs:
        raise ValueError("Nenhum análogo disponível para validação.")

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    params = params or {}
    ranked = sort_by_scientific_score(analogs)
    now = datetime.now().isoformat(timespec="seconds")
    result_hash = _result_hash(analogs)

    files: Dict[str, str] = {}
    files["csv_validado"] = export_validated_csv(ranked, out / "analog_san_resultados_validados.csv")

    files["fig_01_png"] = _make_distribution_png(ranked, out / "fig_01_distribuicoes.png")
    files["fig_01_svg"] = _make_distribution_svg(ranked, out / "fig_01_distribuicoes.svg")
    files["fig_02_png"] = _make_scatter_png(ranked, out / "fig_02_mw_vs_logp.png", "MW", "LogP", "MW vs LogP", "MW (Da)", "LogP")
    files["fig_02_svg"] = _make_basic_svg(ranked, out / "fig_02_mw_vs_logp.svg", "MW vs LogP", "MW", "LogP", "MW (Da)", "LogP")
    files["fig_03_png"] = _make_scatter_png(ranked, out / "fig_03_similaridade_vs_sintese.png", "SimilaridadeTanimoto", "PrioridadeRegra", "Similaridade vs prioridade manual da regra", "Tanimoto", "Prioridade da regra")
    files["fig_03_svg"] = _make_basic_svg(ranked, out / "fig_03_similaridade_vs_sintese.svg", "Similaridade vs prioridade manual da regra", "SimilaridadeTanimoto", "PrioridadeRegra", "Tanimoto", "Prioridade da regra")
    files["fig_04_png"] = _make_radar_png(ranked, out / "fig_04_radar_top5.png")
    files["fig_04_svg"] = _make_radar_svg(ranked, out / "fig_04_radar_top5.svg")
    files["fig_05_png"] = _make_ranking_png(ranked, out / "fig_05_ranking_top10.png")
    files["fig_05_svg"] = _make_ranking_svg(ranked, out / "fig_05_ranking_top10.svg")

    pareto_mode = any("FrentePareto" in analog for analog in ranked)
    ranking_description = (
        "ordem Pareto; ScoreGlobal exibido apenas como métrica auxiliar"
        if pareto_mode else "ordem por ScoreGlobal heurístico"
    )
    captions = [
        "Figura 1. Distribuição de MW, LogP e TPSA dos análogos aprovados.",
        "Figura 2. Dispersão MW vs LogP para avaliação do espaço físico-químico.",
        "Figura 3. Similaridade Tanimoto versus prioridade manual da regra de transformação.",
        f"Figura 4. Perfil multivariado dos cinco candidatos em {ranking_description}.",
        f"Figura 5. Dez candidatos em {ranking_description}.",
    ]
    captions_path = out / "analog_san_legendas_figuras.txt"
    captions_path.write_text("\n".join(captions) + "\n", encoding="utf-8")
    files["legendas"] = str(captions_path)

    report_lines = [
        "ANALOG-SAN - RELATÓRIO DE VALIDAÇÃO DOS ANÁLOGOS",
        "=" * 72,
        f"Data/hora: {now}",
        f"Versão do Analog-San: {APP_VERSION}",
        f"Versão do módulo de validação: {VALIDATION_VERSION}",
        f"SMILES do composto-protótipo: {base_smiles}",
        f"Total de análogos avaliados: {len(ranked)}",
        f"Hash SHA-256 dos resultados: {result_hash}",
        "",
        "Parâmetros de reprodutibilidade:",
        f"- Seed padrão registrada: {params.get('seed', DEFAULT_RANDOM_SEED)}",
        f"- Filtro Veber: {'sim' if params.get('use_veber', True) else 'não'}",
        f"- SA Score máximo: {params.get('max_sa_score', 6.0)} (estimativa heurística)",
        f"- Exploração multietapas: {'sim' if params.get('allow_multistep', False) else 'não'}",
        f"- Rotas: {', '.join(params.get('selected_routes', [])) or 'não informado'}",
        f"- Perfil de objetivos: {json.dumps(params.get('objective_profile', []), ensure_ascii=False) or '[]'}",
        "",
        "Metodologia resumida:",
        "Os análogos foram avaliados por descritores RDKit, regras de Lipinski",
        "e Veber, similaridade Tanimoto por fingerprint Morgan, SA Score,",
        "alertas químicos PAINS/Brenk quando disponíveis e padrões SMARTS",
        "reativos simples. O ScoreGlobal heurístico combina drug-likeness,",
        "faixa físico-química, similaridade, prioridade manual da regra, SA Score",
        "e alertas. A prioridade manual não é estimativa de facilidade de síntese.",
        "As evidências de reação são mostradas por operação; nesta versão, as",
        "transformações ainda não têm precedentes experimentais catalogados.",
        "O SA Score é uma estimativa computacional, não demonstra que exista rota",
        "prática. Planejamento retrosintético, disponibilidade de reagentes,",
        "condições, rendimento, segurança e confirmação experimental não são",
        "avaliados. Alertas PAINS/Brenk são sinalizadores e exigem revisão contextual.",
        "A prioridade da regra não equivale a custo, número de etapas,",
        "rendimento ou facilidade real de síntese. Os aliases históricos",
        "PrioridadeSintese e FacilidadeSintese permanecem só por compatibilidade.",
        "Quando há perfil de objetivos, a ordenação usa frentes de Pareto com",
        "as preferências de descritores informadas, similaridade, SA Score e",
        "alertas; não representa previsão de eficácia ou de síntese.",
        "",
        "Top 10 candidatos (ordem Pareto, se objetivos foram definidos):",
        "ID\tRegra prior.\tClasse hipotética\tTipo registro\tEvidência\tFontes\tRota verificada\tFrente Pareto\tObjetivos\tScore auxiliar\tCategoria\tTanimoto\tSA\tAlertas\tModificação\tSMILES",
    ]
    for analog in ranked[:10]:
        report_lines.append(
            f"{analog.get('ID_Molecula')}\t"
            f"{analog.get('PrioridadeRegra', '')}\t"
            f"{analog.get('EvidenciaClasseReacao', '')}\t"
            f"{analog.get('EvidenciaTipoRegistro', '')}\t"
            f"{analog.get('EvidenciaStatus', '')}\t"
            f"{analog.get('EvidenciaFonteCount', 0)}\t"
            f"{analog.get('EvidenciaRotaVerificada', False)}\t"
            f"{analog.get('FrentePareto', '')}\t"
            f"{analog.get('Objetivos', '')}\t"
            f"{analog.get('ScoreGlobal')}\t"
            f"{analog.get('CategoriaScore')}\t"
            f"{analog.get('SimilaridadeTanimoto')}\t"
            f"{analog.get('SAScore')}\t"
            f"{analog.get('QtdAlertas')}\t"
            f"{analog.get('Modificacao')}\t"
            f"{analog.get('SMILES')}"
        )
    report_lines.extend(["", "Figuras geradas:"])
    report_lines.extend(f"- {os.path.basename(path)}" for key, path in files.items() if key.startswith("fig_"))

    report_path = out / "analog_san_relatorio_validacao.txt"
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    files["relatorio"] = str(report_path)

    metadata = {
        "app": "Analog-San",
        "app_version": APP_VERSION,
        "validation_version": VALIDATION_VERSION,
        "generated_at": now,
        "base_smiles": base_smiles,
        "analog_count": len(ranked),
        "result_sha256": result_hash,
        "params": params,
        "files": files,
    }
    metadata_path = out / "analog_san_metadados_reprodutibilidade.json"
    metadata_path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")
    files["metadados"] = str(metadata_path)

    return {
        "output_dir": str(out),
        "top_id": ranked[0].get("ID_Molecula"),
        "top_score": ranked[0].get("ScoreGlobal"),
        "top_pareto_front": ranked[0].get("FrentePareto"),
        "result_sha256": result_hash,
        "files": files,
    }
