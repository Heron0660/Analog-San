import unittest

import chem_engine as ce


class TierDiverseSelectionTests(unittest.TestCase):
    def test_selects_requested_count_deterministically_and_preserves_tier_order(self):
        kwargs = dict(
            base_smiles="CC(=O)Nc1ccc(O)cc1",
            target_count=5,
            use_veber=True,
            selected_routes=[ce.ROUTE_CLASSIC, ce.ROUTE_NON_CLASSIC],
            random_seed=41242,
            max_sa_score=6.0,
            allow_multistep=False,
            selection_strategy="tier_diverse",
        )
        first = ce.generate_analogs(**kwargs)
        second = ce.generate_analogs(**kwargs)
        self.assertEqual(len(first), 5)
        self.assertEqual([item["SMILES"] for item in first], [item["SMILES"] for item in second])
        priorities = [item["PrioridadeRegra"] for item in first]
        self.assertEqual(priorities, sorted(priorities))
        self.assertEqual([item["OrdemSelecaoDiversa"] for item in first], list(range(1, 6)))

    def test_rejects_conflicting_modes(self):
        with self.assertRaises(ValueError):
            ce.generate_analogs(
                "CCO", 5, allow_multistep=True, selection_strategy="tier_diverse"
            )
        with self.assertRaises(ValueError):
            ce.generate_analogs("CCO", 5, selection_strategy="unknown")


if __name__ == "__main__":
    unittest.main()
