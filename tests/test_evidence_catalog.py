import unittest

from evidence_catalog import validate_catalog


class EvidenceCatalogTests(unittest.TestCase):
    def setUp(self):
        self.registered = {"N-Metilação"}

    def test_empty_catalog_is_valid(self):
        self.assertEqual(validate_catalog({"schema_version": 1, "transformations": {}}, self.registered), {})

    def test_exact_evidence_requires_reaction_and_reviewer(self):
        catalog = {
            "schema_version": 1,
            "transformations": {
                "N-Metilação": {
                    "records": [{
                        "source_id": "doi:10.1000/example",
                        "title": "Example reaction",
                        "url": "https://doi.org/10.1000/example",
                        "scope": "exact_transform",
                    }]
                }
            },
        }
        errors = validate_catalog(catalog, self.registered)
        self.assertIn("reaction_smiles", " ".join(errors["N-Metilação"]))
        self.assertIn("revisor", " ".join(errors["N-Metilação"]))

    def test_rejects_unregistered_names_and_non_https_sources(self):
        catalog = {
            "schema_version": 1,
            "transformations": {
                "N-Metilação": {"records": [{
                    "source_id": "x", "title": "Record", "url": "http://example.org",
                    "scope": "reaction_class",
                }]},
                "invented": {"records": []},
            },
        }
        errors = validate_catalog(catalog, self.registered)
        self.assertIn("HTTPS", " ".join(errors["N-Metilação"]))
        self.assertIn("invented", errors)


if __name__ == "__main__":
    unittest.main()
