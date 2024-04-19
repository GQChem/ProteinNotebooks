#!/usr/bin/bash

echo "A. AlphaFold (localcolabfold)"
if [ ! -d "localcolabfold" ]; then
    echo
    echo "A.i Cloning repository..."
    wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
    echo
    echo "A.ii Installing..."
    bash install_colabbatch_linux.sh
    rm install_colabbatch_linux.sh
else
    echo "Already installed"
fi

echo
echo "B. RFdiffusion"
if [ ! -d "RFdiffusion" ]; then
    echo
    echo "B.i Cloning repository..."
    git clone https://github.com/RosettaCommons/RFdiffusion.git
    echo
    echo "B.ii Downloading models..."
    cd RFdiffusion
    mkdir models && cd models
    wget http://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
    wget http://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
    wget http://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt
    wget http://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt
    wget http://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt
    wget http://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
    wget http://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt
    echo
    echo "B.iii Installing SE3Transformer..."
    cd ../env/SE3Transformer
    pip install --no-cache-dir -r requirements.txt
    python setup.py install
    cd ../..
    pip install -e . 
    cd ..
else
    echo "Already installed"
fi

echo
echo "C. ProteinMPNN"
if [ ! -d "ProteinMPNN" ]; then
    echo
    echo "C.i Cloning repository..."
    git clone https://github.com/dauparas/ProteinMPNN
else
    echo "Already installed"
fi

echo
echo "D. OmegaFold"
if [ ! -d "OmegaFold" ]; then
    git clone https://github.com/HeliXonProtein/OmegaFold
    cd OmegaFold
    python setup.py install
    wget https://helixon.s3.amazonaws.com/release1.pt
    wget https://helixon.s3.amazonaws.com/release2.pt
    cd ..
else
    echo Already installed
fi

echo
echo "Models setup completed!"
