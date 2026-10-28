from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

from evidence_catalog import load_catalog
import chem_engine as ce


def main() -> int:
    catalog_path = ROOT / "evidence_catalog.json"
    names = [name for name, *_ in ce.MODIFICATION_REGISTRY]
    try:
        catalog = load_catalog(str(catalog_path), names)
    except (OSError, ValueError) as exc:
        print(exc, file=sys.stderr)
        return 1
    record_count = sum(len(item["records"]) for item in catalog["transformations"].values())
    print(f"Valid catalogue: {record_count} source records for {len(names)} registered transformations.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
