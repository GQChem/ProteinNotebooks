#!/bin/bash

module load mamba

notebooks=false
models=false

#Pre-process long options and convert them to short options
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

if $notebooks; then
  echo
  echo "Setting up notebooks..."
  if [ -f "protein_notebooks.sh" ]; then
  rm protein_notebooks.sh
  fi
  wget https://raw.githubusercontent.com/GQChem/ProteinNotebooks/main/InstallationFiles/protein_notebooks.sh
  chmod +x protein_notebooks.sh
  bash protein_notebooks.sh
  rm protein_notebooks.sh
fi

# Install ProteinEnv if not present
env_exists=$(conda env list | grep 'ProteinEnv')

if [ -z "$env_exists" ]; then
  cd ProteinNotebooks/InstallationFiles
  mamba env create -f ProteinEnv.yml 
  cd ../..
else
  echo "ProteinEnv found"
fi

conda activate ProteinEnv #Common to all models, contains pymol as well

if [ "$CONDA_DEFAULT_ENV" = "ProteinEnv" ]; then
  echo "Updating environment"
  if ! conda list | grep -q "pymol-bundle"; then
    echo "Installing Pymol..."
    conda install -c conda-forge -c schrodinger pymol-bundle
  else
    echo "Pymol is already installed."
  fi
  if ! pip list | grep -q "prody"; then
    echo "Installing prody..."
    pip install prody
  else
    echo "Prody is already installed"
  fi
  if ! pip list | grep -q "ipython"; then
    echo "Adding ProteinEnv to Jupyter Kernels"
    pip install ipython
    pip install ipykernel
    ipython kernel install --user --name ProteinEnv
  else
    echo "Kernel is already installed"
  fi
  echo "Environment i up-to-date"

  if $models; then
    echo
    echo "Setting up models..."
    cd /home/$USER/data
    if [ -f "protein_models.sh" ]; then
    rm protein_models.sh
    fi
    wget https://raw.githubusercontent.com/GQChem/ProteinNotebooks/main/InstallationFiles/protein_models.sh
    chmod +x protein_models.sh
    bash protein_models.sh
    rm protein_models.sh
  fi

  echo
  echo "Setup completed"

else
  mamba init
  echo "Restart your shell!"
fi
