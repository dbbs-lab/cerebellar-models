#!/bin/bash

cd cerebellum
python mapped_placement.py
status=$?
[ $status -ne 0 ] && echo "\n\nPlease read SETUP.MD. You are missing some steps."
