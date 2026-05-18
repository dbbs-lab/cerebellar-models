import json
import os
import unittest
from collections.abc import Mapping
from os.path import abspath, dirname, join
from unittest.mock import patch

from bsb import parse_configuration_file
from click.testing import CliRunner

from cerebellar_models.cli import (
    CerebOption,
    MicrozonesParams,
    TypeTermElem,
    _filter_simulations_devices,
    _update_cell_types,
    configure,
)


def deep_equal(d, u, path="/"):
    for k, v in u.items():
        if isinstance(v, Mapping):
            deep_equal(d.get(k, {}), v, path + k + "/")
        elif k not in d or d[k] != v:
            return False  # pragma: no cover
    return True


def mock_print_panel(options, title="test"):
    return


class TestCerebOption(unittest.TestCase):
    def test_number_type(self):
        opt = CerebOption("n", "title", type_term=TypeTermElem.Number)
        self.assertEqual(list(opt.choices), [0])

    def test_text_type(self):
        opt = CerebOption("t", "title", type_term=TypeTermElem.Text)
        self.assertEqual(list(opt.choices), [""])

    def test_boolean_type(self):
        opt = CerebOption("b", "title", type_term=TypeTermElem.Boolean)
        self.assertEqual(list(opt.choices), [False, True])

    def test_selection_without_choices_raises(self):
        with self.assertRaises(TypeError):
            CerebOption("s", "title")  # type_term defaults to Selection, choices=None

    def test_basket_without_choices_raises(self):
        with self.assertRaises(TypeError):
            CerebOption("s", "title", type_term=TypeTermElem.Basket)

    def test_selection_default_is_first_choice(self):
        opt = CerebOption("s", "title", choices=["a", "b", "c"])
        self.assertEqual(opt.value, "a")

    def test_basket_default_is_empty_list(self):
        opt = CerebOption("b", "title", choices=["a", "b"], type_term=TypeTermElem.Basket)
        self.assertEqual(list(opt.value), [])

    def test_custom_default_value(self):
        opt = CerebOption("s", "title", choices=["a", "b"], default_value="b")
        self.assertEqual(opt.value, "b")


class TestUpdateCellTypes(unittest.TestCase):
    def test_new_top_level_key(self):
        """A cell type config key absent from the base config is assigned directly."""
        config = {"cell_types": {"granule_cell": {}}}
        cell_type_configs = {"extra": {"new_section": {"key": "value"}}}
        result = _update_cell_types(config, ["extra"], cell_type_configs)
        self.assertIn("new_section", result)
        self.assertEqual(result["new_section"], {"key": "value"})

    def test_deep_updates_existing_key(self):
        config = {"cell_types": {"granule_cell": {}}}
        cell_type_configs = {"extra": {"cell_types": {"dcn_p": {"radius": 9.5}}}}
        result = _update_cell_types(config, ["extra"], cell_type_configs)
        self.assertIn("dcn_p", result["cell_types"])
        self.assertIn("granule_cell", result["cell_types"])

    def test_network_xyz_maxed(self):
        """Network x/y/z values are maxed across multiple cell type configs."""
        config = {"network": {"x": 400, "y": 400, "z": 100}}
        cell_type_configs = {"extra": {"network": {"x": 200, "z": 495}}}
        result = _update_cell_types(config, ["extra"], cell_type_configs)
        self.assertEqual(result["network"]["x"], 400)  # max(400, 200)
        self.assertEqual(result["network"]["z"], 495)  # max(100, 495)


class TestFilterSimulationsDevices(unittest.TestCase):
    @staticmethod
    def _make_simulations(device_name, cell_models):
        return {
            "sim_file": {
                "simulations": {
                    "sim_name": {
                        "devices": {
                            device_name: {
                                "targetting": {
                                    "cell_models": cell_models,
                                    "strategy": "cell_model",
                                }
                            }
                        }
                    }
                }
            }
        }

    def test_filters_simulation_targeting_missing_cell_type(self):
        sims = self._make_simulations("cf_stimulus", ["io"])
        result = _filter_simulations_devices(sims, ["granule_cell"])
        self.assertNotIn("sim_name", result["sim_file"]["simulations"])

    def test_no_filter_when_all_cell_types_present(self):
        sims = self._make_simulations("mf_stimulus", ["mossy_fibers"])
        result = _filter_simulations_devices(sims, ["mossy_fibers", "granule_cell"])
        self.assertIn("sim_name", result["sim_file"]["simulations"])

    def test_record_devices_not_filtered(self):
        """Devices with 'record' in the name are skipped even if cell type is missing."""
        sims = self._make_simulations("io_record", ["io"])
        result = _filter_simulations_devices(sims, ["granule_cell"])
        self.assertIn("sim_name", result["sim_file"]["simulations"])


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
                abspath(__file__),
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
        folder = dirname(__file__)
        config2 = parse_configuration_file(
            join(folder, "test_configurations/canonical_mouse_awake_io_nest.json")
        ).__tree__()
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
        config = parse_configuration_file("./circuit.yaml").__tree__()

        self.assertTrue(deep_equal(config, config2))
        os.remove("./circuit.yaml")

        # Test default parameters
        result = runner.invoke(configure, ["--microzones"])
        self.assertEqual(result.exit_code, 0)
        config = parse_configuration_file("./circuit.yaml").__tree__()

        self.assertTrue(deep_equal(config, config2))
        os.remove("./circuit.yaml")

    @patch(
        "cerebellar_models.cli.print_panel",
        lambda options, title: mock_print_panel(options, title),
    )
    def test_configure_no_microzones(self):
        runner = CliRunner()
        result = runner.invoke(
            configure,
            ["--species", "mouse", "--output_folder", os.getcwd(), "--extension", "yaml"],
        )
        self.assertEqual(result.exit_code, 0)
        config = parse_configuration_file("./circuit.yaml").__tree__()
        self.assertNotIn("after_placement", config)
        os.remove("./circuit.yaml")

    @patch(
        "cerebellar_models.cli.print_panel",
        lambda options, title: mock_print_panel(options, title),
    )
    def test_configure_json_extension(self):
        runner = CliRunner()
        result = runner.invoke(
            configure,
            ["--species", "mouse", "--output_folder", os.getcwd(), "--extension", "json"],
        )
        self.assertEqual(result.exit_code, 0)
        self.assertTrue(os.path.exists("./circuit.json"))
        config = parse_configuration_file("./circuit.json").__tree__()
        self.assertIn("cell_types", config)
        os.remove("./circuit.json")
