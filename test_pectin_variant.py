import re
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
        # chain_index and bead_index are no longer encoded; catalog names are shared
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 3, "PR", 0.9), "PctRep09")
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 10, "PN", 2.1), "PctNeu01")
        self.assertEqual(build_sweep._pectin_atomtype_name(2, 19, "PC", 4.8), "PctXlk09")
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 5, "PC", 4.0), "PctXlk01")
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 1, "PN", 3.5), "PctNeu15")

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
        # Atomtype is now a shared catalog name (no chain/bead position encoding)
        for chain_idx in (1, 2):
            for bead_idx in (1, 30):
                atype = str(assignments[chain_idx][bead_idx]["atomtype"])
                self.assertTrue(
                    atype.startswith("PctRep") or atype.startswith("PctNeu") or atype.startswith("PctXlk"),
                    msg=f"Unexpected atomtype prefix: {atype}",
                )

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
        # Core C/X/P types come first
        self.assertEqual(names[:3], ["C", "X", "P"])
        # All catalog types must be present
        catalog = build_sweep.catalog_type_names()
        self.assertEqual(len(catalog), 51)  # 20 PR + 20 PN + 11 PC
        for cname in catalog:
            self.assertIn(cname, names)

    def test_base_itp_per_bead_epsilon_matches_assignment(self):
        rng = random.Random(99)
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=rng)
        # The assigned atomtype name must be exactly the catalog name for that epsilon
        for a in build_sweep.iter_assignments(assignments):
            expected = build_sweep._catalog_type_name(str(a["bead_type"]), float(a["epsilon"]))
            self.assertEqual(str(a["atomtype"]), expected,
                msg=f"bead_type={a['bead_type']} epsilon={a['epsilon']}")

    def test_base_itp_has_core_nonbond_params(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(1))
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "sudowoodo_base.itp"
            build_sweep.write_per_bead_base_itp(out_path, assignments)
            text = out_path.read_text()
        # The six core C/X/P pairs must be present (file also contains catalog pairs)
        for pair_str in ["C            C", "C            X", "C            P",
                         "X            X", "X            P", "P            P"]:
            self.assertIn(pair_str, text, msg=f"Missing core nonbond pair: {pair_str!r}")

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
                mol_lines = [l for l in itp.splitlines()
                             if re.match(r'^\s+Pctn\s+\d', l)]
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

    def test_per_fiber_itps_use_catalog_atomtypes(self):
        """Every bead across every fiber must use a shared catalog atomtype."""
        assignments = build_sweep.assign_all_chain_bead_epsilons(4, rng=random.Random(77))
        catalog = set(build_sweep.catalog_type_names())
        for a in build_sweep.iter_assignments(assignments):
            self.assertIn(str(a["atomtype"]), catalog,
                msg=f"atomtype {a['atomtype']!r} not in catalog")

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
        """Appending should inject the full catalog into an existing ITP without
        disturbing nonbond_params; a second call must be a no-op."""
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(10))
        with tempfile.TemporaryDirectory() as tmpdir:
            base_path = Path(tmpdir) / "sudowoodo_base.itp"
            # Minimal pre-existing ITP (as afm_build_sweep.py would write)
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
        # Core C/X/P types still present
        self.assertIn("C  1 100.0", text)
        # nonbond_params block still present and unchanged
        self.assertIn("C C 1 2.673 2.5", text)
        # All catalog type names are now in the file
        for cname in build_sweep.catalog_type_names():
            self.assertIn(cname, text, msg=f"catalog type {cname!r} missing")
        # Catalog entries appear BEFORE [ nonbond_params ]
        idx_nb = text.index("[ nonbond_params ]")
        for cname in build_sweep.catalog_type_names():
            self.assertLess(text.index(cname), idx_nb,
                msg=f"{cname!r} not before nonbond_params")
        # Second call must be a no-op (no duplication)
        with tempfile.TemporaryDirectory() as tmpdir2:
            base_path2 = Path(tmpdir2) / "sudowoodo_base.itp"
            base_path2.write_text(text)
            build_sweep.append_per_bead_atomtypes(base_path2, assignments)
            self.assertEqual(base_path2.read_text(), text)

    # ------------------------------------------------------------------ #
    # Catalog naming and cross-term rules
    # ------------------------------------------------------------------ #

    def test_catalog_type_names_count_and_format(self):
        """51 catalog types total: 20 PR + 20 PN + 11 PC, in that order."""
        names = build_sweep.catalog_type_names()
        self.assertEqual(len(names), 51)
        pr_names = [n for n in names if n.startswith("PctRep")]
        pn_names = [n for n in names if n.startswith("PctNeu")]
        pc_names = [n for n in names if n.startswith("PctXlk")]
        self.assertEqual(len(pr_names), 20)
        self.assertEqual(len(pn_names), 20)
        self.assertEqual(len(pc_names), 11)
        # PR comes before PN comes before PC in the returned list
        self.assertEqual(names[:20], pr_names)
        self.assertEqual(names[20:40], pn_names)
        self.assertEqual(names[40:], pc_names)
        # Spot-check specific names
        self.assertIn("PctRep01", names)
        self.assertIn("PctRep20", names)
        self.assertIn("PctNeu01", names)
        self.assertIn("PctNeu20", names)
        self.assertIn("PctXlk01", names)
        self.assertIn("PctXlk11", names)

    def test_catalog_cross_eps_pn_uses_max(self):
        """PN×PN cross-term = max of the two stored offsets (epsilon − 2.0)."""
        # PN.1 (ε=2.1, offset=0.1), PN.2 (ε=2.2, offset=0.2)
        self.assertAlmostEqual(
            build_sweep._catalog_cross_eps("PN", 2.1, "PN", 2.1), 0.1)
        self.assertAlmostEqual(
            build_sweep._catalog_cross_eps("PN", 2.1, "PN", 2.2), 0.2)
        self.assertAlmostEqual(
            build_sweep._catalog_cross_eps("PN", 2.1, "PN", 2.3), 0.3)

    def test_catalog_cross_eps_pc_uses_mean(self):
        """PC×PC cross-term = arithmetic mean of the two absolute epsilons."""
        self.assertAlmostEqual(
            build_sweep._catalog_cross_eps("PC", 4.0, "PC", 4.0), 4.0)
        self.assertAlmostEqual(
            build_sweep._catalog_cross_eps("PC", 4.0, "PC", 4.1), 4.05)
        self.assertAlmostEqual(
            build_sweep._catalog_cross_eps("PC", 4.1, "PC", 4.1), 4.1)

    def test_base_itp_contains_pn_cross_terms(self):
        # PctNeu01–PctNeu01 and PctNeu01–PctNeu02 cross-terms should appear in nonbond_params.
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(1))
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "sudowoodo_base.itp"
            build_sweep.write_per_bead_base_itp(out_path, assignments)
            text = out_path.read_text()
        # PctNeu01 self-interaction: stored ε = 0.1 (offset from 2.0)
        self.assertIn("PctNeu01", text)
        # Check the nonbond_params block contains PctNeu01 PctNeu01 with epsilon 0.1
        nb_start = text.index("[ nonbond_params ]")
        nb_text = text[nb_start:]
        self.assertRegex(nb_text, r"PctNeu01\s+PctNeu01\s+1\s+1\.000000\s+0\.100000")
        self.assertRegex(nb_text, r"PctNeu01\s+PctNeu02\s+1\s+1\.000000\s+0\.200000")

    def test_base_itp_contains_pc_cross_terms(self):
        # PctXlk01–PctXlk01 and PctXlk01–PctXlk02 cross-terms should appear in nonbond_params.
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, rng=random.Random(1))
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "sudowoodo_base.itp"
            build_sweep.write_per_bead_base_itp(out_path, assignments)
            text = out_path.read_text()
        nb_start = text.index("[ nonbond_params ]")
        nb_text = text[nb_start:]
        self.assertRegex(nb_text, r"PctXlk01\s+PctXlk01\s+1\s+1\.000000\s+4\.000000")
        self.assertRegex(nb_text, r"PctXlk01\s+PctXlk02\s+1\s+1\.000000\s+4\.050000")


if __name__ == "__main__":
    unittest.main()
