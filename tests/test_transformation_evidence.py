import unittest
import json
import tempfile
from pathlib import Path

import chem_engine as ce


class TransformationEvidenceTests(unittest.TestCase):
    def test_every_registered_operation_has_an_explicit_evidence_record(self):
        operation_names = {item[0] for item in ce.MODIFICATION_REGISTRY}
        self.assertEqual(operation_names, set(ce.TRANSFORMATION_EVIDENCE))
        for name in operation_names:
            with self.subTest(operation=name):
                evidence = ce.get_transformation_evidence(name)
                self.assertEqual(evidence["status"], "not_catalogued")
                self.assertTrue(evidence["reaction_class"])
                self.assertEqual(evidence["source_records"], ())
                self.assertIsNone(evidence["conditions"])
                self.assertFalse(evidence["route_verified"])

    def test_catalog_import_labels_precedents_without_verifying_routes(self):
        name = ce.MODIFICATION_REGISTRY[0][0]
        catalog = {
            "schema_version": 1,
            "transformations": {
                name: {"records": [{
                    "source_id": "doi:10.1000/example",
                    "title": "Example reaction-class source",
                    "url": "https://doi.org/10.1000/example",
                    "scope": "reaction_class",
                    "conditions": "Reported conditions",
                }]}
            },
        }
        with tempfile.TemporaryDirectory() as temp_dir:
            path = Path(temp_dir) / "catalog.json"
            path.write_text(json.dumps(catalog), encoding="utf-8")
            loaded = ce.load_transformation_evidence_catalog(str(path))[name]
        self.assertEqual(loaded["status"], "reaction_class_precedent")
        self.assertEqual(len(loaded["source_records"]), 1)
        self.assertFalse(loaded["route_verified"])

    def test_generated_candidates_expose_manual_priority_and_evidence_status(self):
        analogs = ce.generate_analogs(
            "CC(=O)Nc1ccc(O)cc1", target_count=3, random_seed=41242,
            allow_multistep=False,
        )
        self.assertTrue(analogs)
        for analog in analogs:
            self.assertEqual(analog["PrioridadeRegra"], analog["PrioridadeSintese"])
            self.assertEqual(analog["NivelRegraManual"], analog["FacilidadeSintese"])
            self.assertEqual(analog["EvidenciaStatus"], "not_catalogued")
            self.assertEqual(analog["EvidenciaFonteCount"], 0)
            self.assertEqual(analog["EvidenciaFontesJSON"], "[]")
            self.assertIn("nao catalogado", analog["EvidenciaObservacao"])


if __name__ == "__main__":
    unittest.main()
