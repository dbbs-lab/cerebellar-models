import copy
import os
from enum import Enum
from os.path import abspath, dirname, join

import click
import numpy as np
import survey
from bsb import Configuration, format_configuration_content, parse_configuration_file

from cerebellum import __version__
from cerebellum.utils import deep_update, get_folders_in_folder, load_configs_in_folder

ROOT_FOLDER = dirname(dirname(abspath(__file__)))
CONFIGURATION_FOLDER = join(ROOT_FOLDER, "configurations")


class TypeTermElem(Enum):
    Selection = 1
    Basket = 2
    Text = 3
    Number = 4
    Boolean = 5


class CerebOption:
    def __init__(self, name, title, choices=None, default_value=None, type=TypeTermElem.Selection):
        self.name = name
        self.title = title
        if type is TypeTermElem.Number:
            choices = [0]
        elif type is TypeTermElem.Text:
            choices = [""]
        elif type is TypeTermElem.Boolean:
            choices = [False, True]
        elif choices is None:
            raise TypeError("Provide a list of choices for Selection or Basket options")
        self.choices = np.array(choices)
        self.value = default_value or (self.choices[0] if type != TypeTermElem.Basket else [])
        self.type = type

    def menu(self):
        if self.type == TypeTermElem.Basket:
            return survey.widgets.Basket(
                options=self.choices, active=np.where(np.isin(self.choices, self.value))[0]
            )
        elif self.type == TypeTermElem.Selection:
            return survey.widgets.Select(options=self.choices)
        elif self.type == TypeTermElem.Text:
            return survey.widgets.Input(value=str(self.value))
        elif self.type == TypeTermElem.Boolean:
            return survey.widgets.Inquire(default=bool(self.value))
        elif self.type == TypeTermElem.Number:
            return survey.widgets.Numeric(value=float(self.value), decimal=True)

    def main_menu_str(self):
        return f"{self.name} [{self.value}]"


def print_panel(options, title="Configure your cerebellum circuit."):
    form = survey.routines.form(title, form={option.name: option.menu() for option in options})
    for option in options:
        if option.type == TypeTermElem.Basket:
            option.value = option.choices[np.array(list(form[option.name]), dtype=int)]
        elif option.type == TypeTermElem.Selection:
            option.value = option.choices[form[option.name]]
        else:
            option.value = form[option.name]


@click.group(help="Cerebellum CLI")
@click.version_option(__version__)
def app():
    """The main CLI entry point"""
    pass


EXISTING_DIR_PATH = click.Path(exists=True, readable=True, dir_okay=True, resolve_path=True)
AVAILABLE_SPECIES = click.Choice(get_folders_in_folder(CONFIGURATION_FOLDER), case_sensitive=True)
AVAILABLE_EXTENSIONS = click.Choice(["yaml", "json"], case_sensitive=True)


def _configure_species(species):
    main_options = [
        CerebOption(
            "Species",
            "Select a species from the following list:",
            AVAILABLE_SPECIES.choices,
            species,
        ),
    ]
    print_panel(main_options, "Select your configuration's species")
    return main_options[0].value


def _configure_cell_types(species_folder, config_cell_types):
    cell_type_names = []
    for filename1, config_1 in config_cell_types.items():
        cell_types1 = list(config_1["cell_types"].keys())
        for filename2, config_2 in config_cell_types.items():
            if filename1 != filename2 and np.all(
                np.isin(cell_types1, list(config_2["cell_types"].keys()))
            ):
                cell_type_names.insert(0, filename1)
                break
        if filename1 not in cell_type_names:
            cell_type_names.append(filename1)

    species_options = [
        CerebOption(
            "State",
            "Select the state of the subject to model from the following list:",
            get_folders_in_folder(species_folder, {"cell_types"}),
            # default_value="awake",
        ),
        CerebOption(
            "Extra cell types",
            "Select the optional cell types that you want in the final configuration from the following list:",
            cell_type_names,
            type=TypeTermElem.Basket,
            # default_value=["dcn", "io"],
        ),
    ]
    print_panel(
        species_options, "Select the state of the subject and the cell types to add in the circuit"
    )
    return [option.value for option in species_options]


def _update_cell_types(configuration, cell_types, config_cell_types):
    for cell_type in cell_types:
        config_ = config_cell_types[cell_type]
        for k, v in config_.items():  # update within the main components
            if k == "network":
                for net_key, net_v in v.items():
                    if net_key in ["x", "y", "z"]:
                        configuration[k][net_key] = max(configuration[k][net_key], net_v)
            else:
                if k not in configuration:
                    configuration[k] = v
                else:
                    deep_update(configuration[k], v)
    return configuration


def _configure_simulations(config_simulations):
    simulation_names = list(
        set(
            [
                f"{simulator}_{simu}"
                for simulator, v in config_simulations.items()
                for file_ in v["simulations"].values()
                for simu in file_["simulations"].keys()
            ]
        )
    )
    simulator_options = [
        CerebOption(
            "Simulations",
            "Select the simulations(s) that you want to perform from the following list:",
            simulation_names,
            type=TypeTermElem.Basket,
            # default_value=["nest_basal_activity"],
        ),
    ]
    print_panel(
        simulator_options, "Select the simulations(s) that you want your circuit to perform"
    )
    return simulator_options[0].value


def _configure_sim_params(config_simulations, simulation_names):
    dict_sim = {"simulations": {}}
    choices = {}
    for sim_name in simulation_names:
        simulator, simulation = sim_name.split("_", 1)
        for k, v in config_simulations[simulator]["simulations"].items():
            if simulation in v["simulations"]:
                for sim, params in v.items():
                    if sim == "simulations":
                        dict_sim[sim][simulation] = params[simulation]
                    else:
                        dict_sim[sim] = params
        simulation_options = [
            CerebOption(
                "Cell models",
                f"Select the model of neuron for the simulation {sim_name} from the following list:",
                list(config_simulations[simulator]["cell_models"].keys()),
                # default_value="eglif_cond_alpha_multisyn"
            ),
            CerebOption(
                "Connection models",
                f"Select the model of synapse for the simulation {sim_name} from the following list:",
                list(config_simulations[simulator]["connection_models"].keys()),
                # default_value="static_synapse",
            ),
        ]
        print_panel(
            simulation_options,
            f"Select the neuron and synapse model to use during the simulation {sim_name}.",
        )
        deep_update(
            dict_sim, config_simulations[simulator]["cell_models"][simulation_options[0].value]
        )
        deep_update(
            dict_sim,
            config_simulations[simulator]["connection_models"][simulation_options[1].value],
        )
        choices[sim_name] = [c.value for c in simulation_options]
        # Add simulator to simulation name so that we avoid duplicates
        dict_sim["simulations"][sim_name] = dict_sim["simulations"][simulation]
        del dict_sim["simulations"][simulation]
    return dict_sim, choices


def _clear_unnecessary_params(configuration):
    dict_cells = {}
    dict_conns = {}
    dict_devices = {}
    for sim_name, simulation in configuration["simulations"].items():
        dict_cells[sim_name] = []
        dict_conns[sim_name] = []
        dict_devices[sim_name] = []
        for cell in simulation["cell_models"]:
            if cell not in configuration["cell_types"]:
                dict_cells[sim_name].append(cell)
        for syn in simulation["connection_models"]:
            found = False
            for strat in configuration["connectivity"]:
                if strat in syn:
                    loc_strat = configuration["connectivity"][strat]
                    simple_conn = (
                        len(loc_strat["presynaptic"]["cell_types"]) == 1
                        and len(loc_strat["postsynaptic"]["cell_types"]) == 1
                    )
                    if simple_conn and strat == syn:
                        found = True
                        break
                    elif simple_conn != (strat == syn):
                        continue
                    cells = syn.split(strat + "_", 1)[1].split("_to_")
                    if (
                        cells[0] in loc_strat["presynaptic"]["cell_types"]
                        and cells[1] in loc_strat["postsynaptic"]["cell_types"]
                    ):
                        found = True
                        break
            if not found:
                dict_conns[sim_name].append(syn)
        for device_name, device in simulation["devices"].items():
            for target in device["targetting"]["cell_models"]:
                if target not in configuration["cell_types"]:
                    dict_devices[sim_name].append(device_name)

    for sim_name, to_remove in dict_cells.items():
        for cell in to_remove:
            del configuration["simulations"][sim_name]["cell_models"][cell]
        for syn in dict_conns[sim_name]:
            del configuration["simulations"][sim_name]["connection_models"][syn]
        for device in dict_devices[sim_name]:
            del configuration["simulations"][sim_name]["devices"][device]

    return configuration


def _write_config(configuration, output_folder, extension):
    output_options = [
        CerebOption(
            "Configuration folder",
            f"Configure the folder in which to put the configuration file",
            default_value=output_folder,
            type=TypeTermElem.Text,
        ),
        CerebOption(
            "File extension",
            "Select an extension from the following list:",
            AVAILABLE_EXTENSIONS.choices,
            extension,
        ),
    ]
    print_panel(
        output_options,
        "Configure the folder in which to put the configuration file and its extension.",
    )
    filename = os.path.join(output_options[0].value, f"circuit.{output_options[1].value}")
    configuration = Configuration.default(**configuration)  # Check that the configuration works
    with open(filename, "w") as outfile:
        outfile.write(format_configuration_content(extension, configuration))
    print(f"Created the BSB configuration file: {filename}")


@app.command(help="Create a BSB configuration file for your cerebellum circuit.")
@click.option(
    "--output_folder",
    type=EXISTING_DIR_PATH,
    required=False,
    default=os.getcwd(),
    help="Path where to write the output configuration file.",
)
@click.option(
    "--species",
    type=AVAILABLE_SPECIES,
    required=False,
    default="mouse",
    help="Species to reconstruct the circuit from.",
)
@click.option(
    "--extension",
    type=AVAILABLE_EXTENSIONS,
    required=False,
    default="yaml",
    help="Extension for the configuration file.",
)
def configure(output_folder: str, species: str, extension: str):
    # Step 1: Species choice
    species = _configure_species(species)
    species_folder = join(CONFIGURATION_FOLDER, species)

    # Step 2: state and cell types choice
    configuration = parse_configuration_file(
        join(species_folder, f"{species}_cerebellar_cortex.yaml")
    ).__tree__()

    config_cell_types = load_configs_in_folder(join(species_folder, "cell_types"))
    state, cell_types = _configure_cell_types(species_folder, config_cell_types)
    configuration = _update_cell_types(configuration, cell_types, config_cell_types)
    state_folder = join(species_folder, state)

    config_simulations = {
        simulator: {
            "cell_models": load_configs_in_folder(join(state_folder, simulator, "cell_models")),
            "connection_models": load_configs_in_folder(
                join(state_folder, simulator, "connection_models")
            ),
            "simulations": load_configs_in_folder(join(state_folder, simulator), recursive=False),
        }
        for simulator in get_folders_in_folder(state_folder)
    }

    # Step 3: Simulation choice
    simulation_names = _configure_simulations(config_simulations)

    # Step 4: Simulation models choice
    dict_sim, sim_choices = _configure_sim_params(config_simulations, simulation_names)
    deep_update(configuration, dict_sim)

    # Step 5: remove unnecessary cells and connections
    configuration = _clear_unnecessary_params(configuration)

    # Step 6: Add stimulus simulation
    if species == "mouse":
        sim_names = list(configuration["simulations"].keys())
        for sim_name in sim_names:
            simulation = configuration["simulations"][sim_name]
            simulator, _ = sim_name.split("_", 1)
            if simulator == "nest":
                default_stim = {
                    "mf_stimulus": {
                        "device": "poisson_generator",
                        "rate": 150,
                        "start": 1200,
                        "stop": 1260,
                        "targetting": {
                            "strategy": "sphere",
                            "radius": 90,
                            "origin": [150, 65, 100],
                            "cell_models": ["mossy_fibers"],
                        },
                        "weight": 1,
                        "delay": 0.1,
                    }
                }
                simulation_name = "nest_mf_stimulus"
                if "io" in cell_types:
                    default_stim["mf_stimulus"] = {
                        "device": "poisson_generator",
                        "rate": 40,
                        "start": 1000,
                        "stop": 1260,
                        "targetting": {"strategy": "cell_model", "cell_models": ["mossy_fibers"]},
                        "weight": 1,
                        "delay": 0.1,
                    }
                    default_stim["cf_stimulus"] = {
                        "device": "poisson_generator",
                        "rate": 500,
                        "start": 1250,
                        "stop": 1260,
                        "targetting": {"strategy": "cell_model", "cell_models": ["io"]},
                        "receptor_type": 1,
                        "weight": 55 if state == "vitro" else 100.0,
                        "delay": 0.1,
                    }
                    simulation_name = "nest_mf_cf_stimulus"
                elif "dcn" in cell_types:
                    default_stim["mf_stimulus"]["targetting"]["origin"][2] = 300

                configuration["simulations"][simulation_name] = copy.deepcopy(
                    configuration["simulations"][sim_name]
                )
                deep_update(configuration["simulations"][simulation_name]["devices"], default_stim)
                sim_choices[simulation_name] = sim_choices[sim_name]
            else:
                raise ValueError(
                    f"Only nest configurations are implemented. Provided simulator: {simulator}"
                )
    else:
        raise ValueError(f"Only mouse configuration are implemented. Provided species: {species}")

    # Step 6: output_folder
    print("\n\nYour choices are:")
    print(f"Species: {species}")
    print(f"State: {state}")
    print(f"Cell types: {cell_types}")
    print("Simulations:")
    for simulation, choices in sim_choices.items():
        print(f"\t{simulation}:")
        print(f"\t\tCell model: {choices[0]}")
        print(f"\t\tSynapse model: {choices[1]}")

    _write_config(configuration, output_folder, extension)
