import json
import os
import sys
import unittest
from collections.abc import Mapping
from os.path import abspath, dirname, join
from unittest.mock import MagicMock, patch

from bsb import parse_configuration_file
from click.testing import CliRunner

from cerebellar_models.cli import (
    CerebOption,
    TypeTermElem,
    _filter_simulations_devices,
    _get_compatible_models,
    _update_cell_types,
    build_nestml,
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
        result = _update_cell_types(config, ["extra"], cell_type_configs, use_atlas=False)
        self.assertIn("new_section", result)
        self.assertEqual(result["new_section"], {"key": "value"})

    def test_deep_updates_existing_key(self):
        config = {"cell_types": {"granule_cell": {}}}
        cell_type_configs = {"extra": {"cell_types": {"dcn_p": {"radius": 9.5}}}}
        result = _update_cell_types(config, ["extra"], cell_type_configs, use_atlas=False)
        self.assertIn("dcn_p", result["cell_types"])
        self.assertIn("granule_cell", result["cell_types"])

    def test_network_xyz_maxed(self):
        """Network x/y/z values are maxed across multiple cell type configs."""
        config = {"network": {"x": 400, "y": 400, "z": 100}}
        cell_type_configs = {"extra": {"network": {"x": 200, "z": 495}}}
        result = _update_cell_types(config, ["extra"], cell_type_configs, use_atlas=False)
        self.assertEqual(result["network"]["x"], 400)  # max(400, 200)
        self.assertEqual(result["network"]["z"], 495)  # max(100, 495)


class TestGetCompatibleModels(unittest.TestCase):
    @staticmethod
    def _make_models(*cell_type_lists):
        """Build a {model_name: config} dict where each model covers the given cell types."""
        return {
            f"model_{i}": {"cell_models": {ct: {} for ct in ct_list}}
            for i, ct_list in enumerate(cell_type_lists)
        }

    def test_all_models_compatible_when_no_extra_cell_types(self):
        models = self._make_models(
            ["granule_cell", "purkinje_cell"],
            ["granule_cell", "purkinje_cell"],
        )
        result = _get_compatible_models(models, ["granule_cell", "purkinje_cell"])
        self.assertCountEqual(result, ["model_0", "model_1"])

    def test_filters_model_missing_extra_cell_type(self):
        models = self._make_models(
            ["granule_cell", "purkinje_cell", "dcn_p"],  # compatible
            ["granule_cell", "purkinje_cell"],  # missing dcn_p → filtered out
        )
        result = _get_compatible_models(models, ["granule_cell", "purkinje_cell", "dcn_p"])
        self.assertEqual(result, ["model_0"])

    def test_cell_types_not_in_any_model_are_ignored(self):
        """Cell types absent from all models (e.g. mossy_fibers as parrot) don't affect filtering."""
        models = self._make_models(
            ["granule_cell"],
            ["granule_cell"],
        )
        # mossy_fibers not in any model → not required, both models remain
        result = _get_compatible_models(models, ["granule_cell", "mossy_fibers"])
        self.assertCountEqual(result, ["model_0", "model_1"])

    def test_empty_circuit_cell_types_returns_all_models(self):
        models = self._make_models(["granule_cell"], ["granule_cell", "dcn_p"])
        result = _get_compatible_models(models, [])
        self.assertCountEqual(result, ["model_0", "model_1"])


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


class TestBuildNestml(unittest.TestCase):
    """Tests for the ``build-nestml`` CLI command."""

    def setUp(self):
        # Inject a mock for build_models so tests don't require NEST/NESTML to be installed.
        self.mock_build = MagicMock()
        mock_module = MagicMock()
        mock_module._build_nest_models = self.mock_build
        sys.modules["cerebellar_models.nest_models.build_models"] = mock_module

    def tearDown(self):
        sys.modules.pop("cerebellar_models.nest_models.build_models", None)

    def test_default_call_passes_redo_true(self):
        """Running without any option calls _build_nest_models(redo=True)."""
        runner = CliRunner()
        result = runner.invoke(build_nestml)
        self.assertEqual(result.exit_code, 0)
        self.mock_build.assert_called_once_with(redo=True)

    def test_build_dir_option_forwarded(self):
        """--build_dir is forwarded as a keyword argument."""
        runner = CliRunner()
        result = runner.invoke(build_nestml, ["--build_dir", "/tmp"])
        self.assertEqual(result.exit_code, 0)
        call_kwargs = self.mock_build.call_args.kwargs
        self.assertTrue(call_kwargs["redo"])
        self.assertEqual(call_kwargs["build_dir"], "/tmp")

    def test_model_dir_option_forwarded(self):
        """--model_dir is resolved to an absolute path and forwarded."""
        runner = CliRunner()
        with runner.isolated_filesystem():
            os.makedirs("models")
            result = runner.invoke(build_nestml, ["--model_dir", "models"])
            self.assertEqual(result.exit_code, 0)
            call_kwargs = self.mock_build.call_args.kwargs
            self.assertTrue(call_kwargs["redo"])
            self.assertTrue(call_kwargs["model_dir"].endswith("models"))

    def test_invalid_model_dir_exits_with_error(self):
        """A non-existent --model_dir causes Click to exit with code 2."""
        runner = CliRunner()
        result = runner.invoke(build_nestml, ["--model_dir", "/nonexistent/path"])
        self.assertEqual(result.exit_code, 2)
        self.assertIn("Error", result.output)
