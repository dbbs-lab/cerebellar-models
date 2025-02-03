import os
from os.path import abspath, basename, dirname, join, splitext

import click
from bsb import Configuration, format_configuration_content, parse_configuration_file
from simple_term_menu import TerminalMenu

from cerebellum import __version__
from cerebellum.utils import deep_update, get_folders_in_folder, load_configs_in_folder

ROOT_FOLDER = dirname(dirname(abspath(__file__)))
CONFIGURATION_FOLDER = join(ROOT_FOLDER, "configurations")


class CerebOption:
    def __init__(self, name, title, choices, default_value=None, is_multi=False, dependencies=None):
        self.name = name
        self.title = title
        self.choices = choices
        self.value = default_value or (choices[0] if not is_multi else [])
        self.is_multi = is_multi
        self.dependencies = dependencies

    def menu(self):
        menu = TerminalMenu(
            self.choices,
            multi_select=self.is_multi,
            show_multi_select_hint=self.is_multi,
            multi_select_select_on_accept=False,
            multi_select_empty_ok=True,
            title=self.title,
        )
        menu.show()
        if self.is_multi:
            choices = menu.chosen_menu_entries or ()
            if self.dependencies is not None:
                for choice in choices:
                    while choice in self.dependencies and self.dependencies[choice] not in choices:
                        choices = (self.dependencies[choice],) + choices
                        choice = self.dependencies[choice]
            self.value = choices
        else:
            self.value = menu.chosen_menu_entry

    def main_menu_str(self):
        return f"{self.name} [{self.value}]"


def print_panel(options):
    while True:
        main_options = [option.main_menu_str() for option in options] + ["", "Confirm"]
        menu = TerminalMenu(
            main_options,
            skip_empty_entries=True,
            title="Configure your cerebellum circuit.\nSelect the option you want to change or confirm:",
        )
        idx = menu.show()
        if idx == len(main_options) - 1:
            break
        else:
            option = options[idx]
            option.menu()


@click.group(help="Cerebellum CLI")
@click.version_option(__version__)
def app():
    """The main CLI entry point"""
    pass


EXISTING_DIR_PATH = click.Path(exists=True, readable=True, dir_okay=True, resolve_path=True)
AVAILABLE_SPECIES = click.Choice(get_folders_in_folder(CONFIGURATION_FOLDER), case_sensitive=True)
AVAILABLE_EXTENSIONS = click.Choice(["yaml", "json"], case_sensitive=True)


def _configure_species(species, extension):
    main_options = [
        CerebOption(
            "Species",
            "Select a species from the following list:",
            AVAILABLE_SPECIES.choices,
            species,
        ),
        CerebOption(
            "File extension",
            "Select an extension from the following list:",
            AVAILABLE_EXTENSIONS.choices,
            extension,
        ),
    ]
    print_panel(main_options)
    return [option.value for option in main_options]


def _configure_cell_types(species_folder, config_cell_types, common_cell_types):
    cell_type_names = []
    dependencies = {}

    for filename, config_ in config_cell_types.items():
        for cell_type in config_["cell_types"]:
            if cell_type not in common_cell_types:
                if "$import" in config_:
                    k = splitext(basename((config_["$import"]["ref"]).split("#")[0]))[0]
                    if k in cell_type_names:
                        cell_type_names.insert(cell_type_names.index(k) + 1, filename)
                    else:
                        cell_type_names.append(filename)
                    dependencies[filename] = k
                else:
                    cell_type_names.insert(0, filename)
                break

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
            is_multi=True,
            dependencies=dependencies,
            # default_value=["dcn", "io"],
        ),
    ]
    print_panel(species_options)
    return [option.value for option in species_options]


def _update_cell_types(configuration, cell_types, config_cell_types):
    for cell_type in cell_types:
        config_ = config_cell_types[cell_type]
        for k, v in config_.items():  # update within the main components
            if k == "network":
                for net_key, net_v in v.items():
                    if net_key in ["x", "y", "z"]:
                        configuration[k][net_key] = max(configuration[k][net_key], net_v)
            elif k == "$import":
                continue
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
            is_multi=True,
            # default_value=["nest_basal_activity"],
        ),
    ]
    print_panel(simulator_options)
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
                f"Select the model of neuron for the simulation {simulation} from the following list:",
                list(config_simulations[simulator]["cell_models"].keys()),
                is_multi=False,
                # default_value="eglif_cond_alpha_multisyn"
            ),
            CerebOption(
                "Connection models",
                f"Select the model of synapse for the simulation {simulation} from the following list:",
                list(config_simulations[simulator]["connection_models"].keys()),
                is_multi=False,
                # default_value="static_synapse",
            ),
        ]
        print_panel(simulation_options)
        deep_update(
            dict_sim, config_simulations[simulator]["cell_models"][simulation_options[0].value]
        )
        deep_update(
            dict_sim,
            config_simulations[simulator]["connection_models"][simulation_options[1].value],
        )
        choices[sim_name] = [c.value for c in simulation_options]
    return dict_sim, choices


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
    species, extension = _configure_species(species, extension)
    species_folder = join(CONFIGURATION_FOLDER, species)

    # Step 2: state and cell types choice
    configuration = parse_configuration_file(
        join(species_folder, f"{species}_cerebellar_cortex.yaml")
    ).__tree__()

    config_cell_types = load_configs_in_folder(join(species_folder, "cell_types"))
    state, cell_types = _configure_cell_types(
        species_folder, config_cell_types, configuration["cell_types"]
    )
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

    # Step 6: output_folder
    print("Your choices are:")
    print(f"Species: {species}")
    print(f"File extension: {extension}")
    print(f"State: {state}")
    print(f"Cell types: {cell_types}")
    print("Simulations:")
    for simulation in simulation_names:
        print(f"\t{simulation}:")
        print(f"\t\tCell model: {sim_choices[simulation][0]}")
        print(f"\t\tSynapse model: {sim_choices[simulation][1]}")

    output_folder = (
        input(f"\nConfigure the folder in which to put the configuration file.\n[{output_folder}]:")
        or output_folder
    )
    filename = os.path.join(output_folder, f"circuit.{extension}")
    configuration = Configuration.default(**configuration)  # Check that the configuration works
    with open(filename, "w") as outfile:
        outfile.write(format_configuration_content(extension, configuration))
    print(f"Created the BSB configuration file: {filename}")


# if __name__ == "__main__":
#     configure(os.getcwd(), "mouse", "yaml")
