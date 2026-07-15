import random
import tempfile
import unittest
from pathlib import Path

import build_sweep

REPO_ROOT = Path(__file__).resolve().parent
PECTIN_TEMPLATE = REPO_ROOT / "toppar_custom" / "sudowoodo_pectin.itp"


class TestPectinVariant(unittest.TestCase):
    def test_assignments_use_tenth_steps(self):
        assignments = build_sweep.assign_all_chain_bead_epsilons(1, 30, rng=random.Random(7))
        for assignment in build_sweep.iter_assignments(assignments):
            epsilon = assignment["epsilon"]
            self.assertEqual(epsilon, round(epsilon, 1))
            self.assertIn(build_sweep.classify_pectin_epsilon(epsilon), assignment["atomtype"])
            self.assertIn(f"e{build_sweep.epsilon_to_step_index(epsilon):02d}", assignment["atomtype"])

    def test_atomtype_names_include_type_and_step(self):
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 3, 0.9), "PRe09c1b3")
        self.assertEqual(build_sweep._pectin_atomtype_name(1, 10, 2.1), "Pe21c1b10")
        self.assertEqual(build_sweep._pectin_atomtype_name(2, 19, 4.8), "PCe48c2b19")

    def test_base_itp_uses_standard_atomtypes_and_nonbond_format(self):
        assignments = {
            1: {
                1: {"chain_index": 1, "bead_index": 1, "epsilon": 3.2, "bead_type": "P", "atomtype": "Pe32c1b1"},
                2: {"chain_index": 1, "bead_index": 2, "epsilon": 0.4, "bead_type": "PR", "atomtype": "PRe04c1b2"},
                3: {"chain_index": 1, "bead_index": 3, "epsilon": 4.7, "bead_type": "PC", "atomtype": "PCe47c1b3"},
            }
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            out_path = Path(tmpdir) / "sudowoodo_base.itp"
            build_sweep.write_per_bead_base_itp(out_path, assignments)
            entries = []
            in_atomtypes = False
            for line in out_path.read_text().splitlines():
                stripped = line.strip()
                if stripped == "[ atomtypes ]":
                    in_atomtypes = True
                    continue
                if in_atomtypes and stripped.startswith("["):
                    break
                if in_atomtypes and stripped and not stripped.startswith(";"):
                    parts = stripped.split()
                    entries.append((parts[0], float(parts[-1])))
            self.assertEqual(entries, [
                ("C", 0.0),
                ("X", 0.0),
                ("P", 0.0),
                ("PN", 0.0),
                ("PR", 0.0),
                ("PC", 0.0),
            ])

            lines = out_path.read_text().splitlines()
            start = lines.index("[ nonbond_params ]") + 1
            nonbond_entries = []
            for line in lines[start:]:
                stripped = line.strip()
                if not stripped or stripped.startswith(";"):
                    continue
                if stripped.startswith("["):
                    break
                nonbond_entries.append(stripped.split())
            self.assertEqual(nonbond_entries, [
                ["C", "C", "1", "2.673000", "1.000000"],
                ["C", "X", "1", "2.087000", "1.000000"],
                ["C", "P", "1", "1.837000", "1.000000"],
                ["X", "X", "1", "1.500000", "1.000000"],
                ["X", "P", "1", "1.250000", "1.000000"],
                ["P", "P", "1", "1.000000", "1.000000"],
            ])

    def test_build_variant_keeps_template_pectin_itp_and_writes_sorted_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            assignments = build_sweep.build_variant(Path(tmpdir), PECTIN_TEMPLATE, rng=random.Random(1))
            pectin_itp = (Path(tmpdir) / "sudowoodo_pectin.itp").read_text()
            template_itp = PECTIN_TEMPLATE.read_text()
            self.assertEqual(pectin_itp, template_itp)
            report_lines = (Path(tmpdir) / "pectin_assignment_report.txt").read_text().splitlines()
            epsilons = [float(line.split()[-1]) for line in report_lines]
            self.assertEqual(epsilons, sorted(epsilons))


if __name__ == "__main__":
    unittest.main()
