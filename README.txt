# smi2cm5
smi2cm5 lets you create atomic descriptors based on CM5 atomic charges computed using semiempirical tight binding (GFN1-xTB).

More information on this method is available in the [RegioML paper](https://doi.org/10.1039/D1DD00032B).

<a href="https://colab.research.google.com/drive/1n3hOlpv2hHXdis66fSq0GAWHjDq0LR0O?usp=sharing">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

-------------------------------------------------------------------------------------------

### INSTALLATION ###

We recommend using anaconda to install the Python 3 environment:
conda env create -f environment.yml && conda activate hecksqm

Then download the binaries of xtb version 6.4.0:
mkdir dep; cd dep; wget https://github.com/grimme-lab/xtb/releases/download/v6.4.0/xtb-210201.tar.xz; tar -xvf ./xtb-210201.tar.xz; cd ..
