import unittest
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import MagicMock, patch

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
            analysis.start_frame_index(0, 0.25)
        with self.assertRaises(ValueError):
            analysis.start_frame_index(10, 0)
        with self.assertRaises(ValueError):
            analysis.start_frame_index(10, 1.5)

    def test_cutoff_angstrom_converts_nm_to_angstrom(self):
        self.assertEqual(analysis.cutoff_angstrom(5.0), 50.0)
        with self.assertRaises(ValueError):
            analysis.cutoff_angstrom(0)

    def test_plot_results_calls_errorbar_with_std_devs(self):
        mock_plt = MagicMock()
        mock_fig = MagicMock()
        mock_plt.figure.return_value = mock_fig

        with TemporaryDirectory() as td:
            root = Path(td)
            with patch.object(analysis, "load_dependencies", return_value=(None, None, mock_plt, None)):
                epsilons = [1.0, 2.0, 3.0]
                averages = [4.0, 5.0, 6.0]
                std_devs = [0.1, 0.2, 0.3]
                analysis.plot_results(root, epsilons, averages, std_devs, 5.0)

        mock_plt.errorbar.assert_called_once_with(
            epsilons, averages, yerr=std_devs, marker="o", capsize=4
        )


if __name__ == "__main__":
    unittest.main()
