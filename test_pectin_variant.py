import random
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import build_sweep as sweep


class PectinVariantTests(unittest.TestCase):
    def test_build_pectin_variant_epsilon_map_uses_exact_requested_pairs_and_defaults_others_to_two(self):
        epsilon_map = sweep.build_pectin_variant_epsilon_map(2.0, -0.5, 5.0)
        self.assertEqual(epsilon_map[("P", "P")], 2.0)
        self.assertEqual(epsilon_map[("P", "PR")], -0.5)
        self.assertEqual(epsilon_map[("P", "PC")], 2.0)
        self.assertEqual(epsilon_map[("PR", "PR")], 2.0)
        self.assertEqual(epsilon_map[("PC", "PC")], 5.0)
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

    def test_scale_epsilon_in_itp_writes_only_pr_and_pc_atomtypes_and_expected_pairs(self):
        with TemporaryDirectory() as td:
            output_path = Path(td) / "scaled.itp"
            sweep.scale_epsilon_in_itp(
                Path(__file__).resolve().parent / "toppar_custom" / "sudowoodo_base.itp",
                output_path,
                sweep.parse_epsilon_map("CC=1.0,CX=1.0,CP=0.7,XX=1.0,XP=0.5,PP=1.0"),
                sweep.build_pectin_variant_epsilon_map(2.0, -0.5, 5.0),
            )
            text = output_path.read_text()
            self.assertNotIn(" PN ", text)
            self.assertIn("C P 1 1.837000 0.700000", text)
            self.assertIn("C PR 1 1.837000 0.700000", text)
            self.assertIn("C PC 1 1.837000 0.700000", text)
            self.assertIn("P P 1 1.000000 1.000000", text)
            self.assertIn("P PR 1 1.000000 -0.500000", text)
            self.assertIn("P PC 1 1.000000 2.000000", text)
            self.assertIn("PR PR 1 1.000000 2.000000", text)
            self.assertIn("PC PC 1 1.000000 5.000000", text)

    def test_get_atom_types_from_itp_reads_p_pr_and_pc_types(self):
        itp = """[atoms]
  1 P 1 Pctn P1 1 0
  2 PR 1 Pctn PR2 2 0
  3 PC 1 Pctn PC3 3 0

[bonds]
"""
        with TemporaryDirectory() as td:
            itp_path = Path(td) / "pectin.itp"
            itp_path.write_text(itp)
            self.assertEqual(
                sweep._get_atom_types_from_itp(itp_path),
                {1: "P", 2: "PR", 3: "PC"},
            )

    def test_write_randomized_pectin_itp_keeps_neutral_p_and_marks_pr_pc_names(self):
        with TemporaryDirectory() as td:
            random.seed(7)
            output_path = Path(td) / "pectin.itp"
            sweep._write_randomized_pectin_itp(
                Path(__file__).resolve().parent / "toppar_custom" / "sudowoodo_pectin.itp",
                output_path,
                "Pctn_1",
            )
            atom_types = sweep._get_atom_types_from_itp(output_path)
            atom_lines = [
                line.split()
                for line in output_path.read_text().splitlines()
                if line.strip() and not line.strip().startswith(';') and line.strip()[0].isdigit()
            ]
            atom_names_by_id = {int(parts[0]): parts[4] for parts in atom_lines[:30]}
            self.assertTrue(any(atom_type == "P" for atom_type in atom_types.values()))
            self.assertTrue(any(atom_type == "PR" for atom_type in atom_types.values()))
            self.assertTrue(any(atom_type == "PC" for atom_type in atom_types.values()))
            for atom_id, atom_type in atom_types.items():
                self.assertEqual(atom_names_by_id[atom_id], f"{atom_type}{atom_id}")

    def test_update_gro_pectin_atomnames_uses_per_fiber_itp_types(self):
        gro = """Test
    6
    1Pctn    P1    1   0.000   0.000   0.000
    1Pctn    P2    2   0.100   0.100   0.100
    1Pctn    P3    3   0.200   0.200   0.200
    2Pctn    P1    1   0.300   0.300   0.300
    2Pctn    P2    2   0.400   0.400   0.400
    2Pctn    P3    3   0.500   0.500   0.500
   1.00000   1.00000   1.00000
"""
        itp_1 = """[atoms]
  1 P 1 Pctn P1 1 0
  2 PR 1 Pctn PR2 2 0
  3 PC 1 Pctn PC3 3 0

[bonds]
"""
        itp_2 = """[atoms]
  1 PC 1 Pctn PC1 1 0
  2 P 1 Pctn P2 2 0
  3 PR 1 Pctn PR3 3 0

[bonds]
"""
        with TemporaryDirectory() as td:
            td = Path(td)
            gro_path = td / "afm_system.gro"
            toppar_dir = td / "toppar_custom"
            toppar_dir.mkdir()
            gro_path.write_text(gro)
            (toppar_dir / "sudowoodo_pectin_1.itp").write_text(itp_1)
            (toppar_dir / "sudowoodo_pectin_2.itp").write_text(itp_2)

            sweep.update_gro_pectin_atomnames(gro_path, toppar_dir)
            sweep.update_gro_pectin_atomnames(gro_path, toppar_dir)

            updated = gro_path.read_text()
            self.assertIn("1Pctn    P1", updated)
            self.assertIn("1Pctn   PR2", updated)
            self.assertIn("1Pctn   PC3", updated)
            self.assertIn("2Pctn   PC1", updated)
            self.assertIn("2Pctn    P2", updated)
            self.assertIn("2Pctn   PR3", updated)

    def test_parse_pectin_atom_index_handles_p_pr_and_pc_names(self):
        self.assertEqual(sweep._parse_pectin_atom_index("P7"), 7)
        self.assertEqual(sweep._parse_pectin_atom_index("PR8"), 8)
        self.assertEqual(sweep._parse_pectin_atom_index("PC9"), 9)
        self.assertIsNone(sweep._parse_pectin_atom_index("X1"))


if __name__ == "__main__":
    unittest.main()
