#!/bin/bash

source /apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/mamba-23.11.0-0-334ztq7i4mzu762ew2x3kbbrrorhe6eg/etc/profile.d/conda.sh
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


if $models; then
  echo
  # Install ProteinEnv if not present
  env_exists=$(conda env list | grep 'ProteinEnv')

  if [ -z "$env_exists" ]; then
    echo "Creating conda environment..."
    cd ProteinNotebooks/InstallationFiles
    mamba env create -f ProteinEnv.yml 
    cd ../..
  else
    echo "ProteinEnv found"
  fi

  conda activate ProteinEnv #Common to all models, contains pymol as well

  if [ "$CONDA_DEFAULT_ENV" = "ProteinEnv" ]; then
    echo "Updating environment ProteinEnv"

    conda_list_output=$(conda list)
    pip_list_output=$(pip list)

    if echo "$conda_list_output" | grep -q "pymol-bundle"; then
      echo "Pymol is already installed"
    else
      echo "Installing Pymol..."
      conda install -c conda-forge -c schrodinger pymol-bundle
    fi
    if python -c "import prody" &> /dev/null; then
      echo "ProDy is already installed"
    else
      echo "ProDy is not installed, installing now..."
      pip install prody
      pip uninstall biopython
      pip install biopython
    fi
    if echo "$pip_list_output" | grep -q "ipython"; then
      echo "Kernel is already installed"
    else
      echo "Adding ProteinEnv to Jupyter Kernels"
      pip install ipython
      pip install ipykernel
      ipython kernel install --user --name ProteinEnv
    fi
    echo "Environment is up-to-date"

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
fi

echo
echo "Setup completed"
