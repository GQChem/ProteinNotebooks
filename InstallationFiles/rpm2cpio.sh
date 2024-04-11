#!/bin/bash

# Use the Python script from the current working directory
PYTHON_SCRIPT="./rpm2cpio.py"

# Pass all arguments to the Python script
python "$PYTHON_SCRIPT" "$@"