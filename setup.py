import setuptools
import os

_findver = "__version__ = "
_rootpath = os.path.join(os.path.dirname(__file__), "bsb_hdf5", "__init__.py")
with open(_rootpath, "r") as f:
    for line in f:
        if _findver in line:
            f = line.find(_findver)
            __version__ = eval(line[line.find(_findver) + len(_findver) :])
            break
    else:
        raise Exception(f"No `__version__` found in '{_rootpath}'.")

with open("README.md", "r") as fh:
    long_description = fh.read()

requires = ["bsb[NEURON,MPI]>=4.0.0a48"]

setuptools.setup(
    name="cerebellum",
    version=__version__,
    author="Robin De Schepper, Alessio Marta, Andrea Di Francescantonio, Egidio D'Angelo, Claudia Casellato",
    author_email="robingilbert.deschepper@unipv.it",
    description="A full scaffold model of the cerebellum, using the BSB.",
    include_package_data=True,
    package_data={"cerebellum": ["*.json"]},
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dbbs-lab/bsb",
    license="GPLv3",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    entry_points={
        "bsb.config.templates": ["cerebellum_templates = cerebellum:templates"],
    },
    python_requires="~=3.8",
    install_requires=requires,
    project_urls={
        "Bug Tracker": "https://github.com/dbbs-lab/cerebellum/issues/",
        "Documentation": "https://cerebellum.readthedocs.io/",
        "Source Code": "https://github.com/dbbs-lab/cerebellum/",
    },
)
