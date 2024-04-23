#!/usr/bin/bash
source /apps/opt/spack/linux-ubuntu20.04-x86_64/gcc-9.3.0/mamba-23.11.0-0-334ztq7i4mzu762ew2x3kbbrrorhe6eg/etc/profile.d/conda.sh
module load mamba

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
    if echo "$pip_list_output" | grep -q "pyppeteer"; then
        echo "pyppeteer is already installed"
    else
        echo "Installing pyppeteer (headless web interface)"
        pip install pyppeteer
    fi
    echo "Environment is up-to-date"
fi