import json
import os
import unittest
from collections.abc import Mapping
from os.path import abspath, dirname, join
from unittest.mock import patch

import yaml
from click.testing import CliRunner

from cerebellar_models.cli import configure


def deep_equal(d, u, path="/"):
    for k, v in u.items():
        if isinstance(v, Mapping):
            deep_equal(d.get(k, {}), v, path + k + "/")
        elif k not in d or d[k] != v:
            return False  # pragma: no cover
    return True


def mock_print_panel(options, title="test"):
    return


class TestCli(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        ROOT_FOLDER = abspath(dirname(dirname(__file__)))
        os.chdir(ROOT_FOLDER)

    @patch(
        "cerebellar_models.cli.print_panel",
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
        "cerebellar_models.cli.print_panel",
        lambda options, title: mock_print_panel(options, title),
    )
    def test_configure(self):
        folder = os.path.dirname(__file__)
        with open(
            os.path.join(folder, "test_configurations/canonical_mouse_awake_io_nest.json"), "r"
        ) as f:
            config2 = json.loads(f.read())
        runner = CliRunner()
        result = runner.invoke(
            configure,
            [
                "--species",
                "mouse",
                "--output_folder",
                os.getcwd(),
                "--extension",
                "yaml",
                "--microzones",
            ],
        )
        self.assertEqual(result.exit_code, 0)
        with open("./circuit.yaml", "r") as f:
            config = yaml.safe_load(f)

        self.assertTrue(deep_equal(config, config2))
        os.remove("./circuit.yaml")

        # Test default parameters
        result = runner.invoke(configure, ["--microzones"])
        self.assertEqual(result.exit_code, 0)
        with open("./circuit.yaml", "r") as f:
            config = yaml.safe_load(f)

        self.assertTrue(deep_equal(config, config2))
        os.remove("./circuit.yaml")
