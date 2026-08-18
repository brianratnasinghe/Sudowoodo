import random
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import build_sweep


class TestPectinVariant(unittest.TestCase):
    def test_catalog_type_names_are_three_fixed_types(self):
        self.assertEqual(build_sweep.catalog_type_names(), ["PctRep", "PctNeu", "PctXlk"])

    def test_atomtype_names_match_bead_class(self):
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 3, "PR", 0.9), "PctRep")
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 10, "PN", 2.1), "PctNeu")
        self.assertEqual(build_sweep._pectin_atomtype_name(2, 19, "PC", 4.8), "PctXlk")

    def test_assignments_use_fixed_type_epsilons(self):
        epsilon_by_type = build_sweep.pectin_epsilon_by_type(0.4, 1.7, 3.9)
        assignments = build_sweep.assign_all_chain_bead_epsilons(
            1,
            rng=random.Random(7),
            epsilon_by_type=epsilon_by_type,
        )
        for assignment in build_sweep.iter_assignments(assignments):
            self.assertEqual(assignment["epsilon"], epsilon_by_type[assignment["bead_type"]])

    def test_fiber_has_fixed_composition(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(3, rng=random.Random(42))
        for chain_index, chain in assignments.items():
            counts = Counter(a["bead_type"] for a in chain.values())
            self.assertEqual(counts["PC"], build_sweep.PC_PER_FIBER, f"chain {chain_index}")
            self.assertEqual(counts["PR"], build_sweep.PR_PER_FIBER, f"chain {chain_index}")
            self.assertEqual(counts["PN"], build_sweep.PN_PER_FIBER, f"chain {chain_index}")
            self.assertEqual(len(chain), build_sweep.BEADS_PER_FIBER, f"chain {chain_index}")

    def test_fiber_bead_order_is_randomised(self):
        a1 = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(1))
        a2 = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(2))
        seq1 = [a1[1][i]["bead_type"] for i in range(1, build_sweep.BEADS_PER_FIBER + 1)]
        seq2 = [a2[1][i]["bead_type"] for i in range(1, build_sweep.BEADS_PER_FIBER + 1)]
        self.assertNotEqual(seq1, seq2)

    def test_base_itp_has_three_catalog_atomtypes(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(1))
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "sudowoodo_base.itp"
            build_sweep.write_per_bead_base_itp(out_path, assignments)
            text = out_path.read_text()
        for atomtype in build_sweep.catalog_type_names():
            self.assertIn(atomtype, text)
        self.assertNotIn("PctRep01", text)
        self.assertNotIn("PctNeu01", text)
        self.assertNotIn("PctXlk01", text)

    def test_base_itp_contains_custom_pectin_interactions(self):
        epsilon_by_type = build_sweep.pectin_epsilon_by_type(0.4, 1.2, 4.8)
        assignments = build_sweep.assign_all_chain_bead_epsilons(
            1,
            rng=random.Random(1),
            epsilon_by_type=epsilon_by_type,
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "sudowoodo_base.itp"
            build_sweep.write_per_bead_base_itp(out_path, assignments, epsilon_by_type=epsilon_by_type)
            text = out_path.read_text()
        self.assertRegex(text, r"PctRep\s+PctRep\s+1\s+1\.000000\s+0\.400000")
        self.assertRegex(text, r"PctNeu\s+PctNeu\s+1\s+1\.000000\s+1\.200000")
        self.assertRegex(text, r"PctXlk\s+PctXlk\s+1\s+1\.000000\s+4\.800000")
        self.assertRegex(text, r"PctRep\s+PctNeu\s+1\s+1\.000000\s+2\.500000")
        self.assertRegex(text, r"PctNeu\s+PctXlk\s+1\s+1\.000000\s+2\.500000")

    def test_per_fiber_itp_molecule_name(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(3, rng=random.Random(5))
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            build_sweep.write_per_fiber_pectin_itps(tmppath, assignments)
            for chain_index in range(1, 4):
                itp = (tmppath / f"sudowoodo_pectin_{chain_index}.itp").read_text()
                self.assertIn(f"Pctn_{chain_index}", itp)

    def test_per_fiber_itp_atom_names_are_bead_types(self):
        """Atom names in [atoms] section must be PR, PN, or PC (not P1..P30)."""
        assignments = build_sweep.assign_all_chain_bead_epsilons(2, rng=random.Random(42))
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            build_sweep.write_per_fiber_pectin_itps(tmppath, assignments)
            for chain_index in range(1, 3):
                itp_text = (tmppath / f"sudowoodo_pectin_{chain_index}.itp").read_text()
                in_atoms = False
                for line in itp_text.splitlines():
                    stripped = line.split(";")[0].strip()
                    if not stripped:
                        continue
                    if stripped.startswith("["):
                        in_atoms = "atoms" in stripped.lower()
                        continue
                    if in_atoms:
                        parts = stripped.split()
                        if len(parts) >= 5:
                            atom_name = parts[4]
                            self.assertIn(atom_name, {"PR", "PN", "PC"},
                                          f"Unexpected atom name {atom_name!r} in chain {chain_index}")

    def test_per_fiber_itps_use_only_shared_catalog_atomtypes(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(4, rng=random.Random(77))
        catalog = set(build_sweep.catalog_type_names())
        for assignment in build_sweep.iter_assignments(assignments):
            self.assertIn(str(assignment["atomtype"]), catalog)

    def test_build_variant_creates_expected_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            build_sweep.build_variant(tmppath, chain_count=2, rng=random.Random(1))
            self.assertTrue((tmppath / "sudowoodo_base.itp").exists())
            self.assertTrue((tmppath / "sudowoodo_pectin_1.itp").exists())
            self.assertTrue((tmppath / "sudowoodo_pectin_2.itp").exists())
            self.assertTrue((tmppath / "pectin_assignment_report.txt").exists())

    def test_assignment_report_preserves_fiber_order(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            build_sweep.build_variant(tmppath, chain_count=1, rng=random.Random(1))
            report_lines = [line for line in (tmppath / "pectin_assignment_report.txt").read_text().splitlines() if line.strip()]
        self.assertEqual(len(report_lines), build_sweep.BEADS_PER_FIBER)
        self.assertTrue(report_lines[0].startswith("chain=1 bead=1 "))
        self.assertTrue(report_lines[-1].startswith(f"chain=1 bead={build_sweep.BEADS_PER_FIBER} "))

    def test_append_per_bead_atomtypes_preserves_core_parameters(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(10))
        with tempfile.TemporaryDirectory() as tmpdir:
            base_path = Path(tmpdir) / "sudowoodo_base.itp"
            base_path.write_text(
                "[ atomtypes ]\n"
                "C  1 100.0 0.000 A 0.0 0.0\n"
                "X  1  50.0 0.000 A 0.0 0.0\n"
                "P  1  26.6 0.000 A 0.0 0.0\n\n"
                "[ nonbond_params ]\n"
                "C C 1 2.673 9.1\n"
                "C X 1 2.087 8.2\n"
                "C P 1 1.837 7.3\n"
                "X X 1 1.500 6.4\n"
                "X P 1 1.250 5.5\n"
                "P P 1 1.000 4.6\n"
                "C PctRep01 1 1.837 3.2\n"
                "X PctRep01 1 1.250 2.3\n"
                "P PctRep01 1 1.000 1.4\n"
            )
            build_sweep.append_per_bead_atomtypes(base_path, assignments)
            text = base_path.read_text()
        self.assertRegex(text, r"C\s+C\s+1\s+2\.673000\s+9\.100000")
        self.assertRegex(text, r"C\s+PctRep\s+1\s+1\.837000\s+3\.200000")
        self.assertRegex(text, r"X\s+PctNeu\s+1\s+1\.250000\s+2\.300000")
        self.assertRegex(text, r"P\s+PctXlk\s+1\s+1\.000000\s+1\.400000")


if __name__ == "__main__":
    unittest.main()
