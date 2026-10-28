import unittest

import validation_engine as ve


class ParetoObjectiveTests(unittest.TestCase):
    def candidate(self, smiles, logp, similarity, sa, alerts, priority, score):
        return {
            "SMILES": smiles,
            "LogP": logp,
            "SimilaridadeTanimoto": similarity,
            "SAScore": sa,
            "QtdAlertas": alerts,
            "PrioridadeSintese": priority,
            "ScoreGlobal": score,
        }

    def test_assigns_nondominated_fronts_and_orders_first_front(self):
        candidates = [
            self.candidate("a", 2.0, 0.9, 2.0, 0, 99, 1),
            self.candidate("b", 3.0, 0.8, 3.0, 0, 1, 99),
            self.candidate("c", 1.0, 0.7, 5.0, 2, 1, 99),
        ]
        ranked = ve.rank_by_objectives(
            candidates, {"LogP": 1.0}, [{"property": "LogP", "direction": "increase"}]
        )
        fronts = {item["SMILES"]: item["FrentePareto"] for item in ranked}
        self.assertEqual(fronts, {"a": 1, "b": 1, "c": 2})
        self.assertEqual(ranked[0]["FrentePareto"], 1)
        self.assertEqual(ranked[-1]["SMILES"], "c")

    def test_maintain_prefers_small_absolute_change(self):
        candidates = [
            self.candidate("near", 2.1, 0.9, 2.0, 0, 1, 99),
            self.candidate("far", 4.0, 0.7, 3.0, 0, 1, 99),
        ]
        ranked = ve.rank_by_objectives(
            candidates, {"LogP": 2.0}, [{"property": "LogP", "direction": "maintain"}]
        )
        fronts = {item["SMILES"]: item["FrentePareto"] for item in ranked}
        self.assertEqual(fronts["near"], 1)
        self.assertGreater(fronts["far"], fronts["near"])
        self.assertAlmostEqual(ranked[0]["DeltaObjetivos"]["LogP"], 0.1)

    def test_rejects_duplicate_and_unknown_objectives(self):
        for profile in (
            [{"property": "LogP", "direction": "increase"}] * 2,
            [{"property": "unknown", "direction": "increase"}],
        ):
            with self.subTest(profile=profile), self.assertRaises(ValueError):
                ve.rank_by_objectives([], {"LogP": 1.0}, profile)


if __name__ == "__main__":
    unittest.main()
