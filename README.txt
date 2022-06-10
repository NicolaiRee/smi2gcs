-------------------------------------------------------------------------------------------

### INSTALLATION ###

We recommend using anaconda to install the Python 3 environment:
conda env create -f environment.yml && conda activate hecksqm

Then download the binaries of xtb version 6.4.0:
mkdir dep; cd dep; wget https://github.com/grimme-lab/xtb/releases/download/v6.4.0/xtb-210201.tar.xz; tar -xvf ./xtb-210201.tar.xz; cd ..
