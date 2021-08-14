#!/bin/sh

# prepare environment
# python3.8 -m venv venv
python -m venv venv
. venv/bin/activate 
cd ..

# get bsb v4.0 outside
git clone https://github.com/dbbs-lab/bsb.git
cd bsb
git checkout v4.0 
pip install -e .[dev]
pre-commit install
cd ..

# get from drive the lingula files
# something like: 
# curl -XGET https://drive.google.com/file/d/11WaZ0cQK0PVI5z1QPNoAQNyBh4mik_VL/view
