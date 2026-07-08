import random
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import afm_build_sweep as sweep


class PectinVariantTests(unittest.TestCase):
    def test_build_pectin_variant_epsilon_map_uses_exact_pn_pairs_and_defaults_others_to_two(self):
        epsilon_map = sweep.build_pectin_variant_epsilon_map(2.0, -0.5, 5.0)
        self.assertEqual(epsilon_map[("PN", "PN")], 2.0)
        self.assertEqual(epsilon_map[("PN", "PR")], -0.5)
        self.assertEqual(epsilon_map[("PN", "PC")], 5.0)
        self.assertEqual(epsilon_map[("PR", "PR")], 2.0)
        self.assertEqual(epsilon_map[("PC", "PC")], 2.0)
        self.assertEqual(epsilon_map[("PR", "PC")], 2.0)

    def test_choose_distributed_positions_requires_at_least_four_beads(self):
        with self.assertRaisesRegex(ValueError, "must have at least 4 beads"):
            sweep._choose_distributed_positions(3)

    def test_choose_distributed_positions_sets_exactly_two_and_two(self):
        random.seed(7)
        neg, pos = sweep._choose_distributed_positions(30)
        self.assertEqual(len(neg), 2)
        self.assertEqual(len(pos), 2)
        self.assertTrue(neg.isdisjoint(pos))
        quarters = {idx // (30 // 4) for idx in neg | pos}
        self.assertEqual(len(quarters), 4)

    def test_count_pectin_fibers_from_gro(self):
        gro = """Test
    4
    1Cell    C1    1   0.000   0.000   0.000
    1Pctn    P1    1   0.000   0.000   0.000
    1Pctn    P2    2   0.100   0.100   0.100
    2Pctn    P1    1   0.200   0.200   0.200
   1.00000   1.00000   1.00000
"""
        with TemporaryDirectory() as td:
            gro_path = Path(td) / "afm_system.gro"
            gro_path.write_text(gro)
            self.assertEqual(sweep.count_pectin_fibers_from_gro(gro_path), 2)

    def test_scale_epsilon_in_itp_writes_new_pectin_variant_types(self):
        with TemporaryDirectory() as td:
            output_path = Path(td) / "scaled.itp"
            sweep.scale_epsilon_in_itp(
                Path(__file__).resolve().parent / "toppar_custom" / "sudowoodo_base.itp",
                output_path,
                sweep.parse_epsilon_map("CC=1.0,CX=1.0,CP=0.7,XX=1.0,XP=0.5,PP=1.0"),
                sweep.build_pectin_variant_epsilon_map(2.0, -0.5, 5.0),
            )
            text = output_path.read_text()
            self.assertIn("C PN 1 1.837000 0.700000", text)
            self.assertIn("C PR 1 1.837000 0.700000", text)
            self.assertIn("C PC 1 1.837000 0.700000", text)
            self.assertIn("PN PN 1 1.000000 2.000000", text)
            self.assertIn("PN PR 1 1.000000 -0.500000", text)
            self.assertIn("PN PC 1 1.000000 5.000000", text)
            self.assertIn("PR PR 1 1.000000 2.000000", text)
            self.assertIn("PC PC 1 1.000000 2.000000", text)


if __name__ == "__main__":
    unittest.main()
