import copy
import os
from enum import Enum
from os.path import abspath, dirname, join

import click
import numpy as np
from bsb import (
    Configuration,
    ConfigurationError,
    format_configuration_content,
    parse_configuration_file,
)

from cerebellar_models import __version__
from cerebellar_models.utils import (
    deep_update,
    get_folders_in_folder,
    load_configs_in_folder,
)

ROOT_FOLDER = dirname(dirname(abspath(__file__)))
CONFIGURATION_FOLDER = join(ROOT_FOLDER, "configurations")


class TypeTermElem(Enum):
    """
    Enum for all supported survey widgets
    """

    Selection = 1
    Basket = 2
    Text = 3
    Number = 4
    Boolean = 5


class CerebOption:
    """
    Class to store the information of an element of a survey form
    """

    def __init__(
        self,
        name: str,
        title: str,
        choices=None,
        default_value=None,
        type_term: TypeTermElem = TypeTermElem.Selection,
    ):
        self.name = name
        """Name of the option"""
        self.title = title
        """Help text for the option"""
        if type_term is TypeTermElem.Number:
            choices = [0]
        elif type_term is TypeTermElem.Text:
            choices = [""]
        elif type_term is TypeTermElem.Boolean:
            choices = [False, True]
        elif choices is None:
            raise TypeError("Provide a list of choices for Selection or Basket options")
        self.choices = np.array(choices)
        """List of possible choices for this option"""
        self.value = default_value or (self.choices[0] if type_term != TypeTermElem.Basket else [])
        """Current value for the option"""
        self.type = type_term
        """Type of the option"""

    def get_widget(self):  # pragma: no cover
        """
        Return the survey widget for this option
        """
        import survey

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


class MicrozonesParams:
    """Class to save the parameters for the micro-zones labelling"""

    def __init__(
        self,
        labels: list[str] = None,
        cell_types: list[str] = None,
    ):
        self.labels = labels or []
        """Labels to associate to each microzone"""
        self.cell_types = cell_types or []
        """Cell types to split into microzones"""
        self.conn_duplicated = {}
        """Duplicated connections obtained from the cell type labelling"""


def print_panel(options, title="Configure your cerebellum circuit."):  # pragma: no cover
    """
    Print a survey form based on a list of options.
    The result of the form will be saved in its options.

    :param list[CerebOption] options: List of options to display
    :param str title: Title to display on top of the form
    """
    import survey

    form = survey.routines.form(
        title, form={option.name: option.get_widget() for option in options}
    )
    for option in options:
        if option.type == TypeTermElem.Basket:
            option.value = option.choices[np.array(list(form[option.name]), dtype=int)]
        elif option.type == TypeTermElem.Selection:
            option.value = option.choices[form[option.name]]
        else:
            option.value = form[option.name]


@click.group(help="Cerebellar models CLI")
@click.version_option(__version__)
def app():  # pragma: no cover
    """The main CLI entry point"""
    pass


EXISTING_DIR_PATH = click.Path(
    exists=True, readable=True, writable=True, file_okay=False, resolve_path=True
)
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


def _configure_cell_types(species_folder, config_cell_types, add_microzones: bool):
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
            default_value="awake",
        ),
        CerebOption(
            "Extra cell types",
            "Select the optional cell types that you want in the final configuration from the following list:",
            cell_type_names,
            type_term=TypeTermElem.Basket,
            default_value=["dcn", "io"],
        ),
        CerebOption(
            "Add microzones?",
            "Split your circuit into microzones",
            type_term=TypeTermElem.Boolean,
            default_value=add_microzones,
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


def _add_microzones(configuration, micro_params: MicrozonesParams):
    # Keep only the labelled cells that are in the final config.
    cells_found = []
    for cell_type in configuration["cell_types"]:
        if cell_type in micro_params.cell_types:
            cells_found.append(cell_type)
    micro_params.cell_types = cells_found

    # add labelling strat
    configuration["after_placement"] = {
        "label_microzones": {
            "strategy": "cerebellar_models.placement.microzones.LabelMicrozones",
            "cell_types": cells_found,
            "labels": micro_params.labels,
        }
    }

    # The connectivity rules that have both their pre- and post- synaptic cell types
    # labelled are updated to filter the first label (no change of name)
    first_label = micro_params.labels[0]
    for strat_name, strat in configuration["connectivity"].items():
        found = True
        for hemitype in ["presynaptic", "postsynaptic"]:
            if not found:
                break
            for cell_type in strat[hemitype]["cell_types"]:
                if cell_type not in micro_params.cell_types:
                    found = False
                    break
        if found:
            duplicate_rules = [copy.deepcopy(strat)] * len(micro_params.labels[1:])
            for hemitype in ["presynaptic", "postsynaptic"]:
                if "labels" not in strat[hemitype]:
                    strat[hemitype]["labels"] = []
                    for duplicate_rule in duplicate_rules:
                        duplicate_rule[hemitype]["labels"] = []

                strat[hemitype]["labels"].append(first_label)
                for label, duplicate_rule in zip(micro_params.labels[1:], duplicate_rules):
                    duplicate_rule[hemitype]["labels"].append(label)
            micro_params.conn_duplicated[strat_name] = duplicate_rules

    # Duplicate the connectivity rule for each extra label
    for strat_name, duplicate_rules in micro_params.conn_duplicated.items():
        for i in range(len(duplicate_rules)):
            new_name = f"{strat_name}_{micro_params.labels[1 + i]}"
            configuration["connectivity"][new_name] = duplicate_rules[i]
            micro_params.conn_duplicated[strat_name][i] = new_name  # keep only new names

    # Update references to the connectivity rules
    for strat_name, strat in configuration["connectivity"].items():
        for k, v in strat.items():
            if isinstance(v, list):
                for elem in v:
                    if elem in micro_params.conn_duplicated:
                        new_names = [f"{elem}_{label}" for label in micro_params.labels[1:]]
                        strat[k].extend(new_names)
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
            type_term=TypeTermElem.Basket,
            default_value=["nest_basal_activity"],
        ),
    ]
    print_panel(
        simulator_options, "Select the simulations(s) that you want your circuit to perform"
    )
    return simulator_options[0].value


def _configure_sim_params(config_simulations, simulation_names, micro_params: MicrozonesParams):
    dict_sim = {"simulations": {}}
    choices = {}
    for sim_name in simulation_names:
        sim_name = str(sim_name)
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
                default_value="tsodyks2_synapse",
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

        # duplicate connections linked to micro-zones
        simulation = dict_sim["simulations"][sim_name]
        for conn, new_names in micro_params.conn_duplicated.items():
            for new_name in new_names:
                simulation["connection_models"][new_name] = copy.deepcopy(
                    simulation["connection_models"][conn]
                )
        # duplicate devices to record labelled cells separately
        names = list(simulation["devices"].keys())
        for cell_type in micro_params.cell_types:
            for name in names:
                device = simulation["devices"][name]
                if (
                    "_record" in name
                    and device["targetting"]["strategy"] == "cell_model"
                    and cell_type in device["targetting"]["cell_models"]
                ):
                    device["targetting"]["strategy"] = "by_label"
                    device["targetting"]["labels"] = [micro_params.labels[0]]
                    for label in micro_params.labels[1:]:
                        new_name = f"{name.split('_record')[0]}_{label}_record"
                        simulation["devices"][new_name] = copy.deepcopy(device)
                        simulation["devices"][new_name]["targetting"]["labels"] = [label]
                    break

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
    output_options = []
    if output_folder is None:
        output_options.append(
            CerebOption(
                "Configuration folder",
                "Configure the folder in which to put the configuration file",
                default_value=os.getcwd(),
                type_term=TypeTermElem.Text,
            )
        )
    else:
        output_folder = abspath(output_folder)
        print(f"Output folder chosen: {output_folder}")
    if extension is None:
        output_options.append(
            CerebOption(
                "File extension",
                "Select an extension from the following list:",
                AVAILABLE_EXTENSIONS.choices,
                default_value=AVAILABLE_EXTENSIONS.choices[0],
            )
        )
    else:
        print(f"File extension chosen: {extension}")
    if len(output_options) > 0:
        print_panel(
            output_options,
            "Configure the folder in which to put the configuration file and its extension.",
        )
        output_folder = output_folder or output_options[0].value
        extension = extension or output_options[-1].value
    filename = os.path.join(output_folder, f"circuit.{extension}")
    try:
        configuration = Configuration.default(**configuration)  # Check that the configuration works
    except ConfigurationError as e:  # pragma: no cover
        raise ValueError("A BSB error happened when loading your cerebellar circuit") from e
    with open(filename, "w") as outfile:
        outfile.write(format_configuration_content(extension, configuration))
    print(f"Created the BSB configuration file: {filename}")


@app.command(help="Create a BSB configuration file for your cerebellum circuit.")
@click.option(
    "--species",
    type=AVAILABLE_SPECIES,
    required=False,
    help="Species to reconstruct the circuit from.",
)
@click.option(
    "--output_folder",
    type=EXISTING_DIR_PATH,
    required=False,
    help="Path where to write the output configuration file.",
)
@click.option(
    "--extension",
    type=AVAILABLE_EXTENSIONS,
    required=False,
    help="Extension for the configuration file.",
)
@click.option(
    "--microzones",
    required=False,
    is_flag=True,
    help="Split your circuit into 2 separated microzones.",
)
def configure(
    species: str = None,
    output_folder: str = None,
    extension: str = None,
    microzones: bool = False,
):
    """
    Resolve a canonical cerebellum configuration file for BSB based on user choices.

    :param str species: species to reconstruct the circuit from.
    :param str output_folder: path where to write the output configuration file.
    :param str extension: extension for the configuration file.
    :param bool microzones: whether to split the circuit into microzones.
    """
    # Step 1: Species choice
    if species is None:
        species = _configure_species(species)
    else:
        print(f"Species chosen: {species}")
    species_folder = join(CONFIGURATION_FOLDER, species)

    # Step 2: state and cell types choice
    configuration = parse_configuration_file(
        join(species_folder, f"{species}_cerebellar_cortex.yaml")
    ).__tree__()

    config_cell_types = load_configs_in_folder(join(species_folder, "cell_types"))
    state, cell_types, microzones = _configure_cell_types(
        species_folder, config_cell_types, microzones
    )
    configuration = _update_cell_types(configuration, cell_types, config_cell_types)
    state_folder = join(species_folder, state)
    if microzones:
        micro_params = MicrozonesParams(
            cell_types=["purkinje_cell", "dcn_p", "dcn_i", "io"], labels=["plus", "minus"]
        )
        configuration = _add_microzones(configuration, micro_params)
    else:
        micro_params = MicrozonesParams()

    # Step 3: Simulation choice
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
    simulation_names = _configure_simulations(config_simulations)

    # Step 4: Simulation models choice
    dict_sim, sim_choices = _configure_sim_params(
        config_simulations, simulation_names, micro_params
    )
    deep_update(configuration, dict_sim)

    # Step 5: remove unnecessary cells and connections
    configuration = _clear_unnecessary_params(configuration)

    # Step 6: Add stimulus simulation
    if species == "mouse":
        sim_names = list(configuration["simulations"].keys())
        for sim_name in sim_names:
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
            else:  # pragma: no cover
                raise ValueError(
                    f"Only nest configurations are implemented. Provided simulator: {simulator}"
                )
    else:  # pragma: no cover
        raise ValueError(f"Only mouse configuration are implemented. Provided species: {species}")

    # Step 7: Recap choices
    print("\n\nYour choices are:")
    print(f"Species: {species}")
    print(f"State: {state}")
    print(f"Cell types: {cell_types}")
    print(f"With microzones: {microzones}")
    if microzones:
        print(f"\tCell types: {micro_params.cell_types}")
        print(f"\tLabels: {micro_params.labels}")
    print("Simulations:")
    for simulation, choices in sim_choices.items():
        print(f"\t{simulation}:")
        print(f"\t\tCell model: {choices[0]}")
        print(f"\t\tSynapse model: {choices[1]}")
    print("----------------------------\n")
    # Step 8: output folder and extension
    _write_config(configuration, output_folder, extension)
