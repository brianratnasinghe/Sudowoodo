import unittest
from argparse import Namespace
from decimal import Decimal
from pathlib import Path
from tempfile import TemporaryDirectory

import build_pectin_monomer_sweep as sweep


class PectinMonomerSweepTests(unittest.TestCase):
    def test_iter_epsilon_values_includes_stop(self):
        self.assertEqual(
            sweep.iter_epsilon_values(Decimal("0.1"), Decimal("0.5"), Decimal("0.1")),
            [
                Decimal("0.1"),
                Decimal("0.2"),
                Decimal("0.3"),
                Decimal("0.4"),
                Decimal("0.5"),
            ],
        )

    def test_generate_positions_returns_requested_count(self):
        box = (40.0, 40.0, 40.0)
        positions = sweep.generate_positions(100, box, 3.5)
        self.assertEqual(len(positions), 100)
        self.assertTrue(
            all(0.0 <= x <= box[0] and 0.0 <= y <= box[1] and 0.0 <= z <= box[2]
                for x, y, z in positions)
        )

    def test_build_case_writes_unbonded_monomer_sweep_files(self):
        args = Namespace(
            count=100,
            box=(40.0, 40.0, 40.0),
            spacing=3.5,
            prod_ns=Decimal("100"),
            dt_ps=Decimal("0.1"),
            gmx="gmx",
            ntomp=24,
            ntmpi=1,
        )
        with TemporaryDirectory() as td:
            case_dir = Path(td) / "pp_eps_0.5"
            sweep.build_case(case_dir, Decimal("0.5"), args)

            gro_text = (case_dir / "afm_system.gro").read_text()
            top_text = (case_dir / "afm_system.top").read_text()
            base_itp = (case_dir / "toppar_custom" / "sudowoodo_base.itp").read_text()
            pectin_itp = (case_dir / "toppar_custom" / "sudowoodo_pectin.itp").read_text()
            production_mdp = (case_dir / "production.mdp").read_text()

            self.assertEqual(gro_text.splitlines()[1].strip(), "100")
            self.assertIn("Pctn 100", top_text)
            self.assertIn("0.500000", base_itp)
            self.assertIn("P     P", base_itp)
            self.assertIn("[atoms]", pectin_itp)
            self.assertNotIn("[bonds]", pectin_itp)
            self.assertIn("nsteps                   = 1000000", production_mdp)


if __name__ == "__main__":
    unittest.main()
