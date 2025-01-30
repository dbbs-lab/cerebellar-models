import json
import os
from collections.abc import Mapping
from glob import glob
from os.path import basename, dirname, isdir, join, splitext

import yaml


def get_folders_in_folder(folder_path: str, excepts: set = None):
    if excepts is None:
        excepts = set()
    return [f for f in os.listdir(folder_path) if isdir(join(folder_path, f)) and f not in excepts]


def load_configs_in_folder(folder_path: str, recursive=True):
    configs = {}

    for ext in ["yaml", "yml", "json"]:
        files = f"/**/*.{ext}" if recursive else f"/*.{ext}"
        configs.update(load_configs_from_files(glob(folder_path + files, recursive=True)))
    return configs


def load_configs_from_files(filenames):
    configs = {}
    loaders = {"yaml": yaml.safe_load, "yml": yaml.safe_load, "json": json.loads}
    for filename in filenames:
        extension = filename.rsplit(".", 1)[1]
        if extension in loaders.keys():
            try:
                with open(filename, "r") as f:
                    data = loaders[extension](f.read())
                    configs[splitext(basename(filename))[0]] = data
            except Exception as e:
                pass
    return configs


def deep_update(d, u):
    for k, v in u.items():
        if isinstance(v, Mapping):
            d[k] = deep_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def find_key(adict, key):
    stack = [[adict, ""]]
    found = []
    while stack:
        d, p = stack.pop()
        if key in d:
            found.append([p, d[key]])
        for k, v in d.items():
            if isinstance(v, dict):
                stack.append([v, p + "/" + k])

    return found


def fetch_reference(content, ref):
    parts = [p for p in ref.split("/")[1:] if p]
    n = content
    loc = ""
    for part in parts:
        loc += "/" + part
        try:
            n = n[part]
        except KeyError:
            raise ValueError(
                "'{}' in File reference '{}' does not exist".format(loc, ref)
            ) from None
        if not isinstance(n, dict):
            raise ValueError(
                "File references can only point to dictionaries. '{}' is a {}".format(
                    "{}' in '{}".format(loc, ref) if loc != ref else ref,
                    type(n).__name__,
                )
            )
    return n


def update_at_ref(orig, content, loc, ref, keys):
    stack_key = [p for p in loc.split("/")[1:] if p]
    while stack_key:
        k = stack_key.pop(0)
        orig = orig[k]
    del orig["$import"]
    c = fetch_reference(content, ref)
    for k in keys:
        if k in orig and isinstance(orig[k], dict):
            old = orig[k]
            deep_update(c[k], old)
        orig[k] = c[k]


def resolve_dependencies(
    configuration_dict,
    folder,
):
    dependencies = {}
    for k, v in configuration_dict.items():
        dependencies[k] = []
        for p, imp in find_key(v, "$import"):
            glob, ref = (imp["ref"]).split("#")
            base = splitext(basename(glob))[0]
            glob = join(folder, glob)
            dependencies[k].append([[glob, base], p, ref, imp["values"]])
    final_dict = {}
    found = []
    for file, deps in dependencies.items():
        for i, dep in enumerate(deps):
            if dep[0][1] == file:
                raise ValueError("Circular dependency detected")
            if dep[0][1] not in configuration_dict.keys():
                deep_update(
                    final_dict,
                    resolve_dependencies(load_configs_from_files([dep[0][0]]), dirname(dep[0][0])),
                )
                found.append(dep[0][1])
            dep[0] = dep[0][1]
            dependencies[file][i] = dep

    while dependencies:
        for file, deps in dependencies.items():
            if not deps:
                deep_update(final_dict, configuration_dict[file])
                found.append(file)
                break
            else:
                to_remove = -1
                for i, dep in enumerate(deps):
                    if dep[0] in found:
                        update_at_ref(configuration_dict[file], final_dict, dep[1], dep[2], dep[3])
                        to_remove = i
                        break
                    elif dep[0] == "" and len(deps) == 1:
                        update_at_ref(configuration_dict[file], final_dict, dep[1], dep[2], dep[3])
                        to_remove = i
                        break
                if to_remove >= 0:
                    del deps[to_remove]
        for f in found:
            if f in dependencies:
                del dependencies[f]
    return final_dict
