import shutil
from os import makedirs
from os.path import abspath, dirname, exists, isdir, join

import appdirs
from bsb.services import MPI
from pynestml.frontend.pynestml_frontend import generate_target

_cereb_dirs = appdirs.AppDirs("cerebellar_models")
_cache_path = _cereb_dirs.user_cache_dir


def _check_nest_is_buildable(module_name):
    """
    Raise a clear error if the installed ``nest`` package cannot be used to build a
    custom NEST extension module (e.g. the NESTML-compiled models in this package).

    ``pynestml``'s builder locates NEST's headers/libraries through the install
    prefix reported by ``nest`` (``nest.build_info["prefix"]``) and the
    ``nest-config`` script found there. The pip-installable ``nest-simulator``
    wheels bundle a self-contained NEST binary without a real install prefix: the
    reported prefix is a leftover build-time temp directory from the wheel's CI
    build, and ``nest-config --includes``/``--libs`` return paths that never
    existed on the local machine. Building against such a NEST install always
    fails, so we detect it upfront instead of surfacing pynestml's cryptic
    ``InvalidPathException``.

    :param str module_name: name of the extension module about to be built, used
        only to make the error message more informative.
    :raises RuntimeError: if the detected NEST install prefix does not exist
        locally, i.e. building a custom extension module is not currently
        possible with this NEST installation.
    """
    import nest

    if "build_info" in dir(nest):
        nest_path = nest.build_info["prefix"]
    else:
        nest_path = nest.ll_api.sli_func("statusdict/prefix ::")
    if not isdir(nest_path):
        raise RuntimeError(
            f"Cannot build the '{module_name}' NEST extension module: the NEST "
            f"install prefix reported by `nest` ('{nest_path}') does not exist on "
            "this machine. This usually happens when NEST was installed from the "
            "'nest-simulator' PyPI wheel (`pip install nest-simulator`), which "
            "bundles a self-contained NEST binary without a real install prefix "
            "and does not currently support building custom NESTML extension "
            "modules (see https://github.com/nest/nest-simulator/issues/3737). "
            "Install NEST from source, conda-forge, or your system package "
            "manager instead, and build the extension module in that environment."
        )


def _build_nest_models(
    model_dir=dirname(__file__),
    build_dir=join(_cache_path, "nest_build"),
    module_name="cerebmodule",
    redo=False,
):
    """
    Build all the nestml models within the provided model directory and deploy them.
    :param str model_dir: Directory containing the nestml files to compile
    :param str build_dir: Directory where the nest models will be compiled
    :param str module_name: Name of the nest module produced as outcome.
    :param bool redo: Flag to force rebuild the nest models.
    """

    model_dir = abspath(model_dir)
    if MPI.get_size() == 1 or MPI.get_rank() == 0:
        if not redo:
            import nest

            nest.ResetKernel()
            try:
                nest.Install(module_name)
                # unload the module
                nest.ResetKernel()
                return
            except nest.NESTErrors.DynamicModuleManagementError as e:
                if "loaded already" in getattr(e, "message", e.args[0]):
                    return
        _check_nest_is_buildable(module_name)
        if not (exists(model_dir) and isdir(model_dir)):
            raise OSError("Model directory does not exist: {}".format(model_dir))
        if exists(build_dir) and isdir(build_dir):
            shutil.rmtree(build_dir)
        makedirs(build_dir)

        generate_target(
            input_path=model_dir,
            target_platform="NEST",
            target_path=build_dir,
            module_name=module_name,
            codegen_opts={
                "gap_junctions": {
                    "enable": True,
                    "membrane_potential_variable": "V_m",
                    "gap_current_port": "I_stim",
                }
            },
        )
