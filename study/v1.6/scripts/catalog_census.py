from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import chem_engine as ce


def main() -> None:
    parser = argparse.ArgumentParser(description="Export the current evidence-catalog census.")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    evidence = ce.TRANSFORMATION_EVIDENCE
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["Transformation", "RegistryRecordType", "SourceRecordCount", "EvidenceStatus", "RouteVerified"])
        for name, *_ in ce.MODIFICATION_REGISTRY:
            item = evidence[name]
            writer.writerow([name, item["record_type"], len(item["source_records"]), item["status"], item["route_verified"]])
    print(f"Wrote census for {len(ce.MODIFICATION_REGISTRY)} transformations to {args.output}")


if __name__ == "__main__":
    main()
