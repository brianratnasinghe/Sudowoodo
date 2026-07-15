import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import analyze_monomer_neighbors as analysis


class AnalyzeMonomerNeighborsTests(unittest.TestCase):
    def test_extract_epsilon_parses_integer_and_decimal_directory_names(self):
        self.assertEqual(analysis.extract_epsilon(Path("pp_eps_4")), 4.0)
        self.assertEqual(analysis.extract_epsilon(Path("pp_eps_0.5")), 0.5)

    def test_iter_case_directories_returns_sorted_pp_eps_directories(self):
        with TemporaryDirectory() as td:
            root = Path(td)
            for name in ("pp_eps_4", "pp_eps_0.5", "pp_eps_2.1", "other"):
                (root / name).mkdir()

            self.assertEqual(
                [path.name for path in analysis.iter_case_directories(root)],
                ["pp_eps_0.5", "pp_eps_2.1", "pp_eps_4"],
            )

    def test_start_frame_index_uses_last_fraction_of_trajectory(self):
        self.assertEqual(analysis.start_frame_index(100, 0.25), 75)
        self.assertEqual(analysis.start_frame_index(5, 0.25), 3)

    def test_start_frame_index_validates_fraction(self):
        with self.assertRaises(ValueError):
            analysis.start_frame_index(10, 0)
        with self.assertRaises(ValueError):
            analysis.start_frame_index(10, 1.5)

    def test_cutoff_angstrom_converts_nm_to_angstrom(self):
        self.assertEqual(analysis.cutoff_angstrom(5.0), 50.0)
        with self.assertRaises(ValueError):
            analysis.cutoff_angstrom(0)


if __name__ == "__main__":
    unittest.main()
