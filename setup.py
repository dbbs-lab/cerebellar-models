import setuptools, sys, os

with open(os.path.join(os.path.dirname(__file__), "cerebellum", "__init__.py"), "r") as f:
    for line in f:
        if "__version__ = " in line:
            exec(line.strip())
            break

with open("README.md", "r") as fh:
    long_description = fh.read()

requires = [
    "bsb[NEURON,MPI]>=4.0.0"
]

setuptools.setup(
    name="cerebellum",
    version=__version__,
    author="Robin De Schepper, Alice Geminiani, Alberto Antonietti, Stefano Casali, Egidio D'Angelo, Claudia Casellato",
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
