#!/bin/bash

# Initialize default values for options
notebooks=false
models=false

# Pre-process long options and convert them to short options
for arg in "$@"; do
  shift
  case "$arg" in
    "--notebooks") set -- "$@" "-n" ;;
    "--models")    set -- "$@" "-m" ;;
    *)            set -- "$@" "$arg"
  esac
done

# Process command-line options
while getopts ":nm" opt; do
  case $opt in
    n)
      notebooks=true
      ;;
    m)
      models=true
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

if $models; then
    echo "Setting up models..."
    echo
    # Check if the file exists
    if [ -f "protein_models.sh" ]; then
    rm protein_models.sh
    fi
    wget https://raw.githubusercontent.com/GQChem/ProteinNotebooks/main/InstallationFiles/protein_models.sh
    chmod +x protein_models.sh
    bash protein_models.sh
    rm protein_models.sh
fi


# Conditional execution based on options
if $notebooks; then
    echo "Setting up notebooks..."
    echo
    # Check if the file exists
    if [ -f "protein_notebooks.sh" ]; then
    rm protein_notebooks.sh
    fi
    wget https://raw.githubusercontent.com/GQChem/ProteinNotebooks/main/InstallationFiles/protein_notebooks.sh
    chmod +x protein_notebooks.sh
    bash protein_notebooks.sh
    rm protein_notebooks.sh
fi

echo
echo Setup completed
