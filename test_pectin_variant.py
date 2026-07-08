import random
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import afm_build_sweep as sweep


class PectinVariantTests(unittest.TestCase):
    def test_choose_distributed_positions_requires_at_least_four_beads(self):
        with self.assertRaises(ValueError):
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
        gro = """Test\n    4\n    1Cell    C1    1   0.000   0.000   0.000\n    1Pctn    P1    1   0.000   0.000   0.000\n    1Pctn    P2    2   0.100   0.100   0.100\n    2Pctn    P1    1   0.200   0.200   0.200\n   1.00000   1.00000   1.00000\n"""
        with TemporaryDirectory() as td:
            gro_path = Path(td) / "afm_system.gro"
            gro_path.write_text(gro)
            self.assertEqual(sweep.count_pectin_fibers_from_gro(gro_path), 2)


if __name__ == "__main__":
    unittest.main()
