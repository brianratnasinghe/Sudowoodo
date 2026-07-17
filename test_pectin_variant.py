import random
import tempfile
import unittest
from collections import Counter
from pathlib import Path

import build_sweep

REPO_ROOT = Path(__file__).resolve().parent


class TestPectinVariant(unittest.TestCase):

    # ------------------------------------------------------------------ #
    # Epsilon-step helpers
    # ------------------------------------------------------------------ #

    def test_epsilon_steps_for_type_pr(self):
        steps = build_sweep._epsilon_steps_for_type("PR")
        self.assertEqual(steps[0], 0.1)
        self.assertEqual(steps[-1], 2.0)
        self.assertEqual(len(steps), 20)

    def test_epsilon_steps_for_type_pn(self):
        steps = build_sweep._epsilon_steps_for_type("PN")
        self.assertEqual(steps[0], 2.1)
        self.assertEqual(steps[-1], 4.0)
        self.assertEqual(len(steps), 20)

    def test_epsilon_steps_for_type_pc(self):
        steps = build_sweep._epsilon_steps_for_type("PC")
        self.assertEqual(steps[0], 4.0)
        self.assertEqual(steps[-1], 5.0)
        self.assertEqual(len(steps), 11)

    def test_classify_pectin_epsilon_boundaries(self):
        self.assertEqual(build_sweep.classify_pectin_epsilon(0.1), "PR")
        self.assertEqual(build_sweep.classify_pectin_epsilon(2.0), "PR")
        self.assertEqual(build_sweep.classify_pectin_epsilon(2.1), "PN")
        self.assertEqual(build_sweep.classify_pectin_epsilon(3.9), "PN")
        # 4.0 is the lower bound of PC and upper bound of PN; PC wins
        self.assertEqual(build_sweep.classify_pectin_epsilon(4.0), "PC")
        self.assertEqual(build_sweep.classify_pectin_epsilon(5.0), "PC")

    # ------------------------------------------------------------------ #
    # Atomtype naming
    # ------------------------------------------------------------------ #

    def test_atomtype_names_include_type_and_step(self):
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 3, "PR", 0.9), "PRe09c1b3")
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 10, "PN", 2.1), "PNe21c1b10")
        self.assertEqual(build_sweep._pectin_atomtype_name(2, 19, "PC", 4.8), "PCe48c2b19")
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 5, "PC", 4.0), "PCe40c1b5")
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 1, "PN", 4.0), "PNe40c1b1")

    # ------------------------------------------------------------------ #
    # Assignment generation
    # ------------------------------------------------------------------ #

    def test_assignments_use_tenth_steps(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(7))
        for assignment in build_sweep.iter_assignments(assignments):
            epsilon = assignment["epsilon"]
            self.assertEqual(epsilon, round(epsilon, 1))
            bead_type = assignment["bead_type"]
            lo, hi = build_sweep.EPSILON_RANGE_BY_TYPE[bead_type]
            self.assertGreaterEqual(epsilon, lo)
            self.assertLessEqual(epsilon, hi)

    def test_fiber_has_fixed_composition(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(3, rng=random.Random(42))
        for chain_index, chain in assignments.items():
            counts = Counter(a["bead_type"] for a in chain.values())
            self.assertEqual(counts["PC"], build_sweep.PC_PER_FIBER, f"chain {chain_index}")
            self.assertEqual(counts["PR"], build_sweep.PR_PER_FIBER, f"chain {chain_index}")
            self.assertEqual(counts["PN"], build_sweep.PN_PER_FIBER, f"chain {chain_index}")
            self.assertEqual(len(chain), build_sweep.BEADS_PER_FIBER, f"chain {chain_index}")

    def test_fiber_bead_order_is_randomised(self):
        # Two chains seeded differently should not have the same bead-type sequence
        a1 = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(1))
        a2 = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(2))
        seq1 = [a1[1][i]["bead_type"] for i in range(1, build_sweep.BEADS_PER_FIBER + 1)]
        seq2 = [a2[1][i]["bead_type"] for i in range(1, build_sweep.BEADS_PER_FIBER + 1)]
        self.assertNotEqual(seq1, seq2)

    def test_assignments_support_multiple_chains(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(2, rng=random.Random(3))
        self.assertEqual(sorted(assignments.keys()), [1, 2])
        self.assertTrue(assignments[1][1]["atomtype"].endswith("c1b1"))
        self.assertTrue(assignments[2][30]["atomtype"].endswith("c2b30"))

    # ------------------------------------------------------------------ #
    # Base ITP generation
    # ------------------------------------------------------------------ #

    def test_base_itp_has_core_atomtypes(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(1))
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "sudowoodo_base.itp"
            build_sweep.write_per_bead_base_itp(out_path, assignments)
            text = out_path.read_text()
            in_atomtypes = False
            names = []
            for line in text.splitlines():
                s = line.strip()
                if s == "[ atomtypes ]":
                    in_atomtypes = True
                    continue
                if in_atomtypes and s.startswith("["):
                    break
                if in_atomtypes and s and not s.startswith(";"):
                    names.append(s.split()[0])
        # Core types must come first
        self.assertEqual(names[:6], ["C", "X", "P", "PN", "PR", "PC"])
        # Per-bead atomtypes follow (30 total for 1 chain)
        per_bead_names = names[6:]
        self.assertEqual(len(per_bead_names), build_sweep.BEADS_PER_FIBER)
        for n in per_bead_names:
            self.assertTrue(
                n.startswith("PR") or n.startswith("PN") or n.startswith("PC"),
                msg=f"Unexpected atomtype prefix: {n}",
            )

    def test_base_itp_per_bead_epsilon_matches_assignment(self):
        rng = random.Random(99)
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=rng)
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "sudowoodo_base.itp"
            build_sweep.write_per_bead_base_itp(out_path, assignments)
            text = out_path.read_text()
            # Build a map atomtype->epsilon from the file
            in_atomtypes = False
            file_map = {}
            for line in text.splitlines():
                s = line.strip()
                if s == "[ atomtypes ]":
                    in_atomtypes = True
                    continue
                if in_atomtypes and s.startswith("["):
                    break
                if in_atomtypes and s and not s.startswith(";"):
                    parts = s.split()
                    if len(parts) >= 7:
                        file_map[parts[0]] = float(parts[-1])
        for a in build_sweep.iter_assignments(assignments):
            atype = str(a["atomtype"])
            self.assertIn(atype, file_map, msg=f"{atype} missing from atomtypes")
            self.assertAlmostEqual(file_map[atype], float(a["epsilon"]), places=5)

    def test_base_itp_has_core_nonbond_params(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(1))
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "sudowoodo_base.itp"
            build_sweep.write_per_bead_base_itp(out_path, assignments)
            lines = out_path.read_text().splitlines()
            start = lines.index("[ nonbond_params ]") + 1
            nonbond_entries = []
            for line in lines[start:]:
                s = line.strip()
                if not s or s.startswith(";"):
                    continue
                if s.startswith("["):
                    break
                parts = s.split()
                if len(parts) == 5:
                    nonbond_entries.append(parts)
        self.assertEqual(nonbond_entries, [
            ["C", "C", "1", "2.673000", "1.000000"],
            ["C", "X", "1", "2.087000", "1.000000"],
            ["C", "P", "1", "1.837000", "1.000000"],
            ["X", "X", "1", "1.500000", "1.000000"],
            ["X", "P", "1", "1.250000", "1.000000"],
            ["P", "P", "1", "1.000000", "1.000000"],
        ])

    # ------------------------------------------------------------------ #
    # Per-fiber ITP generation
    # ------------------------------------------------------------------ #

    def test_per_fiber_itp_molecule_name(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(3, rng=random.Random(5))
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            build_sweep.write_per_fiber_pectin_itps(tmppath, assignments)
            for chain_index in range(1, 4):
                itp = (tmppath / f"sudowoodo_pectin_{chain_index}.itp").read_text()
                self.assertIn(f"Pctn_{chain_index}", itp)
                # The moleculetype name must be Pctn_N, not the bare "Pctn"
                import re as _re
                mol_lines = [l for l in itp.splitlines()
                             if _re.match(r'^\s+Pctn\s+\d', l)]
                self.assertEqual(mol_lines, [],
                    f"Chain {chain_index} has a bare 'Pctn' moleculetype entry")

    def test_per_fiber_itp_atom_count(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(5))
        with tempfile.TemporaryDirectory() as tmpdir:
            build_sweep.write_per_fiber_pectin_itps(Path(tmpdir), assignments)
            itp_text = (Path(tmpdir) / "sudowoodo_pectin_1.itp").read_text()
            in_atoms = False
            atom_count = 0
            for line in itp_text.splitlines():
                s = line.strip()
                if s == "[ atoms ]":
                    in_atoms = True
                    continue
                if in_atoms and s.startswith("["):
                    break
                if in_atoms and s and not s.startswith(";"):
                    atom_count += 1
        self.assertEqual(atom_count, build_sweep.BEADS_PER_FIBER)

    def test_per_fiber_itp_has_bonds_and_angles(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(5))
        with tempfile.TemporaryDirectory() as tmpdir:
            build_sweep.write_per_fiber_pectin_itps(Path(tmpdir), assignments)
            itp_text = (Path(tmpdir) / "sudowoodo_pectin_1.itp").read_text()
        self.assertIn("[ bonds ]", itp_text)
        self.assertIn("[ angles ]", itp_text)
        self.assertIn("k_bond", itp_text)
        self.assertIn("k_theta", itp_text)

    def test_per_fiber_itps_have_distinct_atomtypes(self):
        """Every bead across every fiber must have a unique atomtype name."""
        assignments = build_sweep.assign_all_chain_bead_epsilons(4, rng=random.Random(77))
        all_atomtypes = [a["atomtype"] for a in build_sweep.iter_assignments(assignments)]
        self.assertEqual(len(all_atomtypes), len(set(all_atomtypes)))

    # ------------------------------------------------------------------ #
    # build_variant integration
    # ------------------------------------------------------------------ #

    def test_build_variant_creates_expected_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmppath = Path(tmpdir)
            build_sweep.build_variant(tmppath, chain_count=2, rng=random.Random(1))
            self.assertTrue((tmppath / "sudowoodo_base.itp").exists())
            self.assertTrue((tmppath / "sudowoodo_pectin_1.itp").exists())
            self.assertTrue((tmppath / "sudowoodo_pectin_2.itp").exists())
            self.assertFalse((tmppath / "sudowoodo_pectin.itp").exists(),
                             "Generic single-fiber ITP should not be created")
            self.assertTrue((tmppath / "pectin_assignment_report.txt").exists())

    def test_build_variant_writes_sorted_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            build_sweep.build_variant(Path(tmpdir), rng=random.Random(1))
            report_lines = [
                line for line in (Path(tmpdir) / "pectin_assignment_report.txt").read_text().splitlines()
                if line.strip() and not line.lstrip().startswith(";")
            ]
            epsilons = [float(line.split()[-1]) for line in report_lines]
            self.assertEqual(epsilons, sorted(epsilons))

    # ------------------------------------------------------------------ #
    # append_per_bead_atomtypes
    # ------------------------------------------------------------------ #

    def test_append_per_bead_atomtypes_extends_existing_itp(self):
        """Appending should add per-bead lines without touching nonbond_params."""
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(10))
        with tempfile.TemporaryDirectory() as tmpdir:
            base_path = Path(tmpdir) / "sudowoodo_base.itp"
            # Write a minimal pre-existing base ITP (as afm_build_sweep.py would)
            base_path.write_text(
                "[ atomtypes ]\n"
                "C  1 100.0 0.000 A 0.0 0.0\n"
                "X  1  50.0 0.000 A 0.0 0.0\n"
                "P  1  26.6 0.000 A 0.0 0.0\n"
                "[ nonbond_params ]\n"
                "C C 1 2.673 2.5\n"
            )
            build_sweep.append_per_bead_atomtypes(base_path, assignments)
            text = base_path.read_text()
        # Core types still present
        self.assertIn("C  1 100.0", text)
        # nonbond_params still present and unchanged
        self.assertIn("C C 1 2.673 2.5", text)
        # Per-bead atomtypes present
        for a in build_sweep.iter_assignments(assignments):
            self.assertIn(str(a["atomtype"]), text)
        # Per-bead lines inserted BEFORE [ nonbond_params ]
        idx_nb = text.index("[ nonbond_params ]")
        for a in build_sweep.iter_assignments(assignments):
            self.assertLess(text.index(str(a["atomtype"])), idx_nb)


if __name__ == "__main__":
    unittest.main()
