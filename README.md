# SMI2GCS
With SMI2GCS you can generate atomic descriptors from SMILES.
The atomic descriptors are based on convolutions of CM5 atomic charges computed using semiempirical tight binding (GFN1-xTB).

<p align="center">
  <img src="image/GCS.png" height="400"/>
</p>

More information about the method is available in the [RegioML paper](https://doi.org/10.1039/D1DD00032B).
Including: [1](https://doi.org/10.1002/cmdc.201700097), [2](https://doi.org/10.1002/cmdc.201800309), [3](https://doi.org/10.1002/minf.201800115), [4](https://doi.org/10.1186/s13321-019-0381-4), and [5](https://doi.org/10.1021/acs.jcim.8b00758).

<a href="https://colab.research.google.com/drive/1n3hOlpv2hHXdis66fSq0GAWHjDq0LR0O?usp=sharing">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>


## Installation

We recommend using anaconda to install the Python 3 environment:

    conda env create -f environment.yml && conda activate smi2gcs

Then download the binaries of xtb version 6.4.0:

    mkdir dep; cd dep; wget https://github.com/grimme-lab/xtb/releases/download/v6.4.0/xtb-210201.tar.xz; tar -xvf ./xtb-210201.tar.xz; cd ..


## Details on the sorting algorithm

Sort each shell according to a modified version of the Cahn-Ingold-Prelog (CIP) priority rules and the CM5 charges if CIP is unambiguous:

1) Sort according to atomic number in descending order.
2) If (1) is not unique, for each atom with the same priority (A*):
    <ol><li>Go to bound and yet not included atoms and sum up atomic numbers. Set the priority of A* according to summed atomic numbers.</li>
    <li>If (2i) did not give an unambiguous result expand the shell of each atom A* by one bond.</li>
    <li>Repeat (2ii) until a unique order is found.</li></ol>
3) If no unique order is found in (2) and all bound atoms are included, then sort atoms according to the CM5 charges in descending order.


## Citation 
```
@article{Ree2022,
  title = {RegioML: predicting the regioselectivity of electrophilic aromatic substitution reactions using machine learning},
  volume = {1},
  ISSN = {2635-098X},
  url = {http://dx.doi.org/10.1039/D1DD00032B},
  DOI = {10.1039/d1dd00032b},
  number = {2},
  journal = {Digital Discovery},
  publisher = {Royal Society of Chemistry (RSC)},
  author = {Nicolai Ree and Andreas H. G\"{o}ller and Jan H. Jensen},
  year = {2022},
  pages = {108–114}
}
```
