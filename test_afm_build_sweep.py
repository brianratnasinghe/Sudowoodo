import unittest
from pathlib import Path
from unittest.mock import patch

import afm_build_sweep


class AfmBuildSweepTests(unittest.TestCase):
    @patch("afm_build_sweep.subprocess.run")
    def test_build_afm_system_forwards_pr_pc_bead_counts(self, mock_run):
        afm_build_sweep.build_afm_system(
            seed=42,
            out_dir=Path("."),
            pr_epsilon=0.5,
            pn_epsilon=2.2,
            pc_epsilon=4.5,
            pr_per_fiber=10,
            pc_per_fiber=10,
        )

        cmd = mock_run.call_args.kwargs["args"] if "args" in mock_run.call_args.kwargs else mock_run.call_args.args[0]
        self.assertIn("--pr", cmd)
        self.assertIn("10", cmd)
        self.assertIn("--pc", cmd)
        self.assertEqual(cmd[cmd.index("--pr") + 1], "10")
        self.assertEqual(cmd[cmd.index("--pc") + 1], "10")


if __name__ == "__main__":
    unittest.main()
