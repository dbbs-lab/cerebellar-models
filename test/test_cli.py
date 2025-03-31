import os
import unittest
from collections.abc import Mapping
from unittest.mock import patch

import yaml
from click.testing import CliRunner

from cerebellum.cli import configure


def deep_equal(d, u, path="/"):
    for k, v in u.items():
        if isinstance(v, Mapping):
            deep_equal(d.get(k, {}), v, path + k + "/")
        elif k not in d or d[k] != v:
            return False
    return True


def mock_print_panel(options, title="test"):
    return


class TestCli(unittest.TestCase):
    @patch(
        "cerebellum.cli.print_panel",
        lambda options, title: mock_print_panel(options, title),
    )
    def test_configure_errors(self):
        runner = CliRunner()
        # wrong species
        result = runner.invoke(
            configure, ["--species", "bla", "--output_folder", os.getcwd(), "--extension", "yaml"]
        )
        self.assertEqual(result.exit_code, 2)
        self.assertIn("Error: Invalid value for '--species'", result.output)

        # folder with no write permission
        result = runner.invoke(
            configure, ["--species", "mouse", "--output_folder", "/", "--extension", "yaml"]
        )
        self.assertEqual(result.exit_code, 2)
        self.assertIn("Error: Invalid value for '--output_folder'", result.output)

        # file instead of folder
        result = runner.invoke(
            configure,
            [
                "--species",
                "mouse",
                "--output_folder",
                os.path.abspath(__file__),
                "--extension",
                "yaml",
            ],
        )
        self.assertEqual(result.exit_code, 2)
        self.assertIn("Error: Invalid value for '--output_folder'", result.output)

        # wrong extension
        result = runner.invoke(
            configure, ["--species", "mouse", "--output_folder", os.getcwd(), "--extension", "test"]
        )
        self.assertEqual(result.exit_code, 2)
        self.assertIn("Error: Invalid value for '--extension'", result.output)

    @patch(
        "cerebellum.cli.print_panel",
        lambda options, title: mock_print_panel(options, title),
    )
    def test_configure(self):
        with open("./test/test_configurations/canonical_mouse_awake_io_nest.yaml", "r") as f:
            config2 = yaml.safe_load(f)
        runner = CliRunner()
        result = runner.invoke(
            configure, ["--species", "mouse", "--output_folder", os.getcwd(), "--extension", "yaml"]
        )
        self.assertEqual(result.exit_code, 0)
        with open("./circuit.yaml", "r") as f:
            config = yaml.safe_load(f)

        self.assertTrue(deep_equal(config, config2))
        os.remove("./circuit.yaml")

        # Test default parameters
        result = runner.invoke(configure, [])
        self.assertEqual(result.exit_code, 0)
        with open("./circuit.yaml", "r") as f:
            config = yaml.safe_load(f)

        self.assertTrue(deep_equal(config, config2))
        os.remove("./circuit.yaml")
