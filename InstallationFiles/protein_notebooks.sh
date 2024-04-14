#!/usr/bin/bash

echo
if [ ! -d "ProteinNotebooks" ]; then
    git clone https://github.com/GQChem/ProteinNotebooks
else
    if [ -d "_BackupPDBs" ]; then
        rm -r _BackupPDBs
    fi
    if [ -d "_BackupPipelines" ]; then
        rm -r _BackupPipelines
    fi
    mv ProteinNotebooks/PDBs _BackupPDBs
    mv ProteinNotebooks/Pipelines _BackupPipelines
    rm -rf ProteinNotebooks
    git clone https://github.com/GQChem/ProteinNotebooks
    rm -rf ProteinNotebooks/PDBs
    rm -rf ProteinNotebooks/Pipelines
    mv _BackupPDBs ProteinNotebooks/PDBs
    mv _BackupPipelines ProteinNotebooks/Pipelines
    rm -rf _BackupPDBs
    rm -rf _BackupPipeline
fi

echo
echo "Notebooks setup completed!"
