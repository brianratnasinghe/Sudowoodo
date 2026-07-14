import random
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import build_sweep as sweep


class PectinVariantTests(unittest.TestCase):
    # ------------------------------------------------------------------
    # _choose_distributed_positions (unchanged)
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # count_pectin_fibers_from_gro (unchanged)
    # ------------------------------------------------------------------

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

    # ------------------------------------------------------------------
    # _extract_base_pectin_type
    # ------------------------------------------------------------------

    def test_extract_base_pectin_type_strips_chain_bead_suffix(self):
        self.assertEqual(sweep._extract_base_pectin_type("P1b5"),   "P")
        self.assertEqual(sweep._extract_base_pectin_type("PR2b7"),  "PR")
        self.assertEqual(sweep._extract_base_pectin_type("PC3b15"), "PC")
        # Legacy bare names round-trip unchanged
        self.assertEqual(sweep._extract_base_pectin_type("P"),  "P")
        self.assertEqual(sweep._extract_base_pectin_type("PR"), "PR")
        self.assertEqual(sweep._extract_base_pectin_type("PC"), "PC")

    # ------------------------------------------------------------------
    # _pectin_atomtype_name
    # ------------------------------------------------------------------

    def test_pectin_atomtype_name_encodes_chain_and_bead(self):
        self.assertEqual(sweep._pectin_atomtype_name(1, 5, sweep.PECTIN_NEUTRAL_TYPE),   "P1b5")
        self.assertEqual(sweep._pectin_atomtype_name(2, 7, sweep.PECTIN_REPULSIVE_TYPE), "PR2b7")
        self.assertEqual(sweep._pectin_atomtype_name(3, 15, sweep.PECTIN_CROSSLINK_TYPE), "PC3b15")
        # Large chain/bead indices
        self.assertEqual(sweep._pectin_atomtype_name(2750, 30, sweep.PECTIN_NEUTRAL_TYPE), "P2750b30")

    # ------------------------------------------------------------------
    # assign_all_chain_bead_epsilons
    # ------------------------------------------------------------------

    def test_assign_all_chain_bead_epsilons_type_counts_and_ranges(self):
        random.seed(42)
        chain_beads = sweep.assign_all_chain_bead_epsilons(pectin_count=5, beads_per_chain=30)
        self.assertEqual(len(chain_beads), 5)
        for bead_map in chain_beads:
            self.assertEqual(len(bead_map), 30)
            types = [t for t, _ in bead_map.values()]
            self.assertEqual(types.count(sweep.PECTIN_REPULSIVE_TYPE), 2)
            self.assertEqual(types.count(sweep.PECTIN_CROSSLINK_TYPE), 2)
            neutral_count = types.count(sweep.PECTIN_NEUTRAL_TYPE)
            self.assertEqual(neutral_count, 26)
            for bead_type, eps in bead_map.values():
                lo, hi = sweep.EPSILON_RANGE_BY_TYPE[bead_type]
                self.assertGreaterEqual(eps, lo)
                self.assertLessEqual(eps, hi)

    def test_assign_all_chain_bead_epsilons_each_bead_unique_epsilon(self):
        random.seed(99)
        chain_beads = sweep.assign_all_chain_bead_epsilons(pectin_count=3, beads_per_chain=30)
        all_eps = [eps for bead_map in chain_beads for _t, eps in bead_map.values()]
        # Extremely unlikely that any two randomly drawn floats are identical
        self.assertEqual(len(all_eps), len(set(all_eps)))

    # ------------------------------------------------------------------
    # write_per_bead_base_itp
    # ------------------------------------------------------------------

    def test_write_per_bead_base_itp_structure_and_atomtypes(self):
        random.seed(11)
        chain_beads = sweep.assign_all_chain_bead_epsilons(pectin_count=2, beads_per_chain=30)
        epsilon_map = sweep.parse_epsilon_map("CC=1.2,CX=0.9,CP=1.0,XX=0.8,XP=0.7,PP=1.0")
        with TemporaryDirectory() as td:
            dst = Path(td) / "sudowoodo_base.itp"
            sweep.write_per_bead_base_itp(dst, chain_beads, epsilon_map)
            text = dst.read_text()

            # Structural sections are present
            self.assertIn("[ defaults ]", text)
            self.assertIn("[ atomtypes ]", text)
            self.assertIn("[ nonbond_params ]", text)

            # C and X backbone types are present with real sigmas
            self.assertIn(f"{'C'}", text)
            self.assertIn(f"{'X'}", text)

            # Explicit nonbond_params for C-C, C-X, X-X with user values
            self.assertIn("C      C      1", text)
            self.assertIn("C      X      1", text)
            self.assertIn("X      X      1", text)
            self.assertIn("1.200000", text)   # epsilon_cc = 1.2
            self.assertIn("0.900000", text)   # epsilon_cx = 0.9
            self.assertIn("0.800000", text)   # epsilon_xx = 0.8

            # Every chain×bead atomtype is present
            for chain_idx, bead_map in enumerate(chain_beads, start=1):
                for bead_id, (bead_type, eps) in bead_map.items():
                    name = sweep._pectin_atomtype_name(chain_idx, bead_id, bead_type)
                    self.assertIn(name, text)
                    self.assertIn(f"{eps:.6f}", text)

            # No old monolithic P/PR/PC single-epsilon lines in atomtypes
            lines = text.splitlines()
            atomtype_lines = []
            in_at = False
            for ln in lines:
                if ln.strip().startswith("[ atomtypes ]"):
                    in_at = True
                    continue
                if in_at and ln.strip().startswith("["):
                    break
                if in_at:
                    atomtype_lines.append(ln)
            atomtype_names = []
            for ln in atomtype_lines:
                stripped = ln.strip()
                if stripped and not stripped.startswith(";"):
                    atomtype_names.append(stripped.split()[0])
            self.assertNotIn("P",  atomtype_names)
            self.assertNotIn("PR", atomtype_names)
            self.assertNotIn("PC", atomtype_names)

    # ------------------------------------------------------------------
    # _get_atom_types_from_itp – returns base types even with unique names
    # ------------------------------------------------------------------

    def test_get_atom_types_from_itp_reads_base_types_from_unique_names(self):
        itp = """[atoms]
  1 P1b1 1 Pctn P1 1 0
  2 PR1b7 1 Pctn PR7 2 0
  3 PC1b15 1 Pctn PC15 3 0

[bonds]
"""
        with TemporaryDirectory() as td:
            itp_path = Path(td) / "pectin.itp"
            itp_path.write_text(itp)
            self.assertEqual(
                sweep._get_atom_types_from_itp(itp_path),
                {1: "P", 2: "PR", 3: "PC"},
            )

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

    # ------------------------------------------------------------------
    # _write_randomized_pectin_itp – unique atomtypes per bead
    # ------------------------------------------------------------------

    def test_write_randomized_pectin_itp_uses_unique_per_bead_atomtypes(self):
        random.seed(7)
        chain_beads = sweep.assign_all_chain_bead_epsilons(pectin_count=1, beads_per_chain=30)
        bead_assignments = chain_beads[0]
        with TemporaryDirectory() as td:
            output_path = Path(td) / "pectin.itp"
            sweep._write_randomized_pectin_itp(
                Path(__file__).resolve().parent / "toppar_custom" / "sudowoodo_pectin.itp",
                output_path,
                "Pctn_1",
                1,  # chain_idx
                bead_assignments,
            )
            atom_types = sweep._get_atom_types_from_itp(output_path)
            text = output_path.read_text()

            # Base types are P, PR, PC
            self.assertTrue(any(t == "P"  for t in atom_types.values()))
            self.assertTrue(any(t == "PR" for t in atom_types.values()))
            self.assertTrue(any(t == "PC" for t in atom_types.values()))

            # Exactly 2 repulsive and 2 crosslink per chain
            self.assertEqual(sum(1 for t in atom_types.values() if t == "PR"), 2)
            self.assertEqual(sum(1 for t in atom_types.values() if t == "PC"), 2)

            # Every atom's TYPE field in the ITP is the globally unique name
            for bead_id, bead_type in atom_types.items():
                expected_atomtype = sweep._pectin_atomtype_name(1, bead_id, bead_type)
                self.assertIn(expected_atomtype, text)

            # Atom NAME field (col 4) is the short form P{i}, PR{i}, PC{i}
            atom_lines = [
                ln.split()
                for ln in text.splitlines()
                if ln.strip() and not ln.strip().startswith(";") and ln.strip()[0].isdigit()
            ]
            atom_names_by_id = {int(p[0]): p[4] for p in atom_lines[:30]}
            for atom_id, base_type in atom_types.items():
                self.assertEqual(atom_names_by_id[atom_id], f"{base_type}{atom_id}")

    # ------------------------------------------------------------------
    # _parse_pectin_atom_index (unchanged)
    # ------------------------------------------------------------------

    def test_parse_pectin_atom_index_handles_p_pr_and_pc_names(self):
        self.assertEqual(sweep._parse_pectin_atom_index("P7"), 7)
        self.assertEqual(sweep._parse_pectin_atom_index("PR8"), 8)
        self.assertEqual(sweep._parse_pectin_atom_index("PC9"), 9)
        self.assertIsNone(sweep._parse_pectin_atom_index("X1"))

    # ------------------------------------------------------------------
    # update_gro_pectin_atomnames – still works with unique atomtype names
    # ------------------------------------------------------------------

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
        # ITP files use unique atomtype names in the TYPE column but short
        # names in the ATOM NAME column – _get_atom_types_from_itp extracts
        # the base type so the GRO update still produces P1/PR2/PC3 names.
        itp_1 = """[atoms]
  1 P1b1 1 Pctn P1 1 0
  2 PR1b2 1 Pctn PR2 2 0
  3 PC1b3 1 Pctn PC3 3 0

[bonds]
"""
        itp_2 = """[atoms]
  1 PC2b1 1 Pctn PC1 1 0
  2 P2b2 1 Pctn P2 2 0
  3 PR2b3 1 Pctn PR3 3 0

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


if __name__ == "__main__":
    unittest.main()
