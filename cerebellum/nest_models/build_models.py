import shutil
from os import makedirs
from os.path import exists, isdir, join, realpath

import appdirs
from pynestml.frontend.pynestml_frontend import generate_target

_cereb_dirs = appdirs.AppDirs("cerebellum")
_cache_path = _cereb_dirs.user_cache_dir


def _build_nest_models(
    model_dir=realpath("./"), build_dir=join(_cache_path, "nest_build"), module_name="cerebmodule"
):
    """
    Build all the nestml models within the provided model directory and deploy them.

    :param str model_dir: Directory containing the nestml files to compile
    :param str build_dir: Directory where the nest models will be compiled
    :param str module_name: Name of the nest module produced as outcome.
    """
    if not (exists(model_dir) and isdir(model_dir)):
        raise OSError("Model directory does not exist: {}".format(model_dir))
    if exists(build_dir) and isdir(build_dir):
        shutil.rmtree(build_dir)
    makedirs(build_dir)

    generate_target(
        input_path=model_dir, target_platform="NEST", target_path=build_dir, module_name=module_name
    )


if __name__ == "__main__":
    _build_nest_models()
