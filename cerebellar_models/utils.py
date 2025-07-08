import os
from collections.abc import Mapping
from glob import glob
from os.path import basename, isdir, join, splitext

from bsb import parse_configuration_content_to_dict


def get_folders_in_folder(folder_path: str, excepts: set = None):
    """
    List all folders in a given folder path.

    :param str folder_path: folder path
    :param set[str] excepts: sub folders to ignore
    :return: list of folders
    :rtype: list[str]
    """
    if excepts is None:
        excepts = set()
    return [f for f in os.listdir(folder_path) if isdir(join(folder_path, f)) and f not in excepts]


def load_configs_in_folder(folder_path: str, recursive=True):
    """
    Load all BSB configuration as dictionaries in a given folder path.

    :param str folder_path: folder path
    :param bool recursive: if True will load also configurations in sub folders to ignore
    :return: dict linking each filename to its BSB configuration dictionary
    :rtype: dict
    """
    configs = {}

    for ext in ["yaml", "yml", "json"]:
        files = f"/**/*.{ext}" if recursive else f"/*.{ext}"
        configs.update(load_configs_from_files(glob(folder_path + files, recursive=True)))
    return configs


def load_configs_from_files(filenames):
    """
    Load all BSB configuration as dictionaries for a list of given filenames.

    :param list[str] filenames: list of configuration files
    :return: dict linking each filename to its BSB configuration dictionary
    :rtype: dict
    """
    configs = {}
    for filename in filenames:
        with open(filename, "r") as f:
            data = f.read()
            configs[splitext(basename(filename))[0]] = parse_configuration_content_to_dict(
                data, path=filename
            )[0]
    return configs


def deep_update(d: dict, u: Mapping):
    """
    Recursively update a dictionary d based on a dictionary u.
    u will overwrite any keys also in d.

    :return: merged dictionary
    :rtype: dict
    """
    for k, v in u.items():
        if isinstance(v, Mapping):
            d[k] = deep_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d
