"""Validation and loading for manually curated transformation evidence."""

import json
from pathlib import Path
from typing import Dict, Iterable
from urllib.parse import urlparse


ALLOWED_SCOPES = {"exact_transform", "reaction_class"}
REQUIRED_SOURCE_FIELDS = {"source_id", "title", "url", "scope"}


def validate_catalog(data: Dict, registered_transformations: Iterable[str]) -> Dict[str, list[str]]:
    """Return validation errors grouped by transformation name."""
    errors: Dict[str, list[str]] = {}
    if not isinstance(data, dict) or data.get("schema_version") != 1:
        return {"catalog": ["schema_version deve ser 1 e o catálogo deve ser um objeto JSON."]}

    transformations = data.get("transformations")
    if not isinstance(transformations, dict):
        return {"catalog": ["'transformations' deve ser um objeto JSON."]}

    registered = set(registered_transformations)
    for name in transformations.keys() - registered:
        errors[name] = ["Transformação não existe no registro do Analog-San."]

    for name, entry in transformations.items():
        if name not in registered:
            continue
        issues = []
        if not isinstance(entry, dict) or not isinstance(entry.get("records"), list):
            errors[name] = ["A entrada deve conter uma lista 'records'."]
            continue
        source_ids = set()
        for index, record in enumerate(entry["records"]):
            prefix = f"records[{index}]"
            if not isinstance(record, dict):
                issues.append(f"{prefix} deve ser um objeto.")
                continue
            missing = REQUIRED_SOURCE_FIELDS - record.keys()
            if missing:
                issues.append(f"{prefix} sem campos: {', '.join(sorted(missing))}.")
                continue
            source_id = record["source_id"]
            if not isinstance(source_id, str) or not source_id.strip():
                issues.append(f"{prefix}.source_id não pode ser vazio.")
            elif source_id in source_ids:
                issues.append(f"{prefix}.source_id duplicado: {source_id}.")
            source_ids.add(source_id)
            url = record["url"]
            parsed = urlparse(url) if isinstance(url, str) else None
            if not parsed or parsed.scheme != "https" or not parsed.netloc:
                issues.append(f"{prefix}.url deve ser um endereço HTTPS válido.")
            if not isinstance(record["title"], str) or not record["title"].strip():
                issues.append(f"{prefix}.title não pode ser vazio.")
            if record["scope"] not in ALLOWED_SCOPES:
                issues.append(f"{prefix}.scope deve ser exact_transform ou reaction_class.")
            if record["scope"] == "exact_transform":
                reaction_smiles = record.get("reaction_smiles")
                if not isinstance(reaction_smiles, str) or reaction_smiles.count(">") != 2:
                    issues.append(f"{prefix}: evidência exata exige reaction_smiles no formato reação.")
                if not record.get("reviewer"):
                    issues.append(f"{prefix}: evidência exata exige identificação de revisor.")
            for optional_text in ("conditions", "yield", "substrate_scope", "limitations", "reviewer"):
                value = record.get(optional_text)
                if value is not None and not isinstance(value, str):
                    issues.append(f"{prefix}.{optional_text} deve ser texto.")
        if issues:
            errors[name] = issues
    return errors


def load_catalog(path: str, registered_transformations: Iterable[str]) -> Dict:
    """Load a UTF-8 JSON catalog, rejecting malformed or unsupported records."""
    catalog_path = Path(path)
    with catalog_path.open("r", encoding="utf-8") as stream:
        data = json.load(stream)
    errors = validate_catalog(data, registered_transformations)
    if errors:
        details = "; ".join(f"{name}: {' '.join(messages)}" for name, messages in errors.items())
        raise ValueError(f"Catálogo de evidências inválido: {details}")
    return data
