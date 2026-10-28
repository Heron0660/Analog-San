import unittest

from rdkit import Chem

import chem_engine


class ArylHalideSubstitutionTests(unittest.TestCase):
    parent_smiles = "COc1cc(CCN)c(OC)cc1Br"

    def setUp(self):
        self.parent = Chem.MolFromSmiles(self.parent_smiles)

    def test_alkyl_substitutions_replace_aromatic_bromine(self):
        cases = (
            (chem_engine.replace_aryl_halide_with_methyl, 1),
            (chem_engine.replace_aryl_halide_with_ethyl, 2),
            (chem_engine.replace_aryl_halide_with_isopropyl, 3),
            (chem_engine.replace_aryl_halide_with_tert_butyl, 4),
        )
        parent_carbon_count = sum(a.GetAtomicNum() == 6 for a in self.parent.GetAtoms())
        for transform, added_carbons in cases:
            with self.subTest(transform=transform.__name__):
                products = transform(self.parent)
                self.assertEqual(len(products), 1)
                product = products[0]
                self.assertIsNotNone(Chem.MolFromSmiles(Chem.MolToSmiles(product)))
                self.assertFalse(any(a.GetAtomicNum() == 35 for a in product.GetAtoms()))
                self.assertEqual(
                    sum(a.GetAtomicNum() == 6 for a in product.GetAtoms()),
                    parent_carbon_count + added_carbons,
                )
                self.assertEqual(
                    sum(a.GetAtomicNum() == 8 for a in product.GetAtoms()),
                    sum(a.GetAtomicNum() == 8 for a in self.parent.GetAtoms()),
                )

    def test_fluorinated_aryl_is_not_treated_as_coupling_handle(self):
        fluorinated = Chem.MolFromSmiles("Fc1ccccc1")
        for transform in (
            chem_engine.replace_aryl_halide_with_methyl,
            chem_engine.replace_aryl_halide_with_ethyl,
            chem_engine.replace_aryl_halide_with_isopropyl,
            chem_engine.replace_aryl_halide_with_tert_butyl,
        ):
            self.assertEqual(transform(fluorinated), [])

    def test_tert_butyl_analog_is_available_in_direct_generation(self):
        analogs = chem_engine.generate_analogs(
            self.parent_smiles,
            target_count=100,
            random_seed=41242,
            max_sa_score=6.0,
            allow_multistep=False,
        )
        expected = Chem.MolToSmiles(
            chem_engine.replace_aryl_halide_with_tert_butyl(self.parent)[0]
        )
        matching = [analog for analog in analogs if analog["SMILES"] == expected]
        self.assertEqual(len(matching), 1)
        self.assertEqual(matching[0]["Modificacao"], "Ar-Haleto->terc-Butil")
        self.assertEqual(matching[0]["EtapaGeracao"], "protótipo")

    def test_objective_mode_selects_from_complete_direct_candidate_set(self):
        profile = [{"property": "LogP", "direction": "increase"}]
        all_candidates = chem_engine.generate_analogs(
            self.parent_smiles, target_count=100, random_seed=41242,
            max_sa_score=6.0, allow_multistep=False,
        )
        expected = chem_engine.ve.rank_by_objectives(
            [dict(item) for item in all_candidates],
            chem_engine.calculate_properties(self.parent), profile,
        )[0]["SMILES"]
        selected = chem_engine.generate_analogs(
            self.parent_smiles, target_count=1, random_seed=41242,
            max_sa_score=6.0, allow_multistep=False, objective_profile=profile,
        )
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0]["SMILES"], expected)
        self.assertIn("FrentePareto", selected[0])


if __name__ == "__main__":
    unittest.main()
