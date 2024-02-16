# MIT License
#
# Copyright (c) 2024 Nicolai Ree
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem

from DescriptorCreator.NodeDescGenerator import NodeDescGenerator

# xTB path and calc setup
XTBHOME = '/Users/nicolairee/opt/miniforge3/envs/alfabet'
XTBPATH = os.path.join(XTBHOME, 'share/xtb')
MANPATH = os.path.join(XTBHOME, 'share/man')
LD_LIBRARY_PATH = os.path.join(XTBHOME, 'lib')

# path = os.getcwd()
# XTBHOME = os.path.join(path, 'dep/xtb-6.4.0')
# XTBPATH = os.path.join(path, 'dep/xtb-6.4.0/share/xtb')
# MANPATH = os.path.join(path, 'dep/xtb-6.4.0/share/man')
# LD_LIBRARY_PATH = os.path.join(path, 'dep/xtb-6.4.0/lib')

OMP_NUM_THREADS = '1'
MKL_NUM_THREADS = '1'


class GraphChargeShell():
    """
    Class to generate atomic descriptors from SMILES.
    """

    def __init__(self):
       
        # Make seperate directory for descritptor calculations
        self.SQMroot = self._make_SQMroot()
        if not os.path.exists(self.SQMroot):
            os.mkdir(self.SQMroot)

        # Set env parameters for xTB        
        global XTBHOME
        global XTBPATH
        global MANPATH
        global LD_LIBRARY_PATH
        global OMP_NUM_THREADS
        global MKL_NUM_THREADS
        os.environ['XTBHOME'] = XTBHOME
        os.environ['XTBPATH'] = XTBPATH
        os.environ['MANPATH'] = MANPATH
        os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH
        os.environ["OMP_NUM_THREADS"] = OMP_NUM_THREADS
        os.environ['MKL_NUM_THREADS'] = MKL_NUM_THREADS


    def _make_SQMroot(self):
        """
        Make a pathname for the SQM calculations (GFN1-xTB)
        :return: SQMroot
        """
        cwd = os.getcwd()
        SQMroot = cwd + '/' + 'calculations'
        print(f'SQM folder is: \n{SQMroot}')
        return SQMroot


    def generate_3Dxyz(self, smi, name):

        # Smiles to RDKit mol object
        self.rdkit_mol = Chem.AddHs(Chem.MolFromSmiles(smi))

        # Define new mol
        new_mol = Chem.AddHs(Chem.MolFromSmiles(smi))

        # Set conf parameters
        rot_bond = Chem.rdMolDescriptors.CalcNumRotatableBonds(self.rdkit_mol)
        n_conformers = min(1 + 3 * rot_bond, 20)

        # Embed mol object to get cartesian coordinates
        ps = AllChem.ETKDGv3()
        ps.randomSeed = 123
        ps.useSmallRingTorsions=True
        try:
            embed_out = AllChem.EmbedMultipleConfs(self.rdkit_mol, numConfs=n_conformers, params=ps)
        except:
            embed_out = -1
        if embed_out == -1 or self.rdkit_mol.GetNumConformers() == 0:
            print(f'1st embed failed for {name} with SMILES: {smi}; will try useRandomCoords=True')
            ps = AllChem.ETKDGv3()
            ps.randomSeed = 123
            ps.useSmallRingTorsions=True
            ps.useRandomCoords=True #added 21/6 - 2021
            # ps.maxIterations=1000 #added 21/6 - 2021
            try:
                embed_out = AllChem.EmbedMultipleConfs(self.rdkit_mol, numConfs=n_conformers, params=ps)
            except:
                embed_out = -1
            if embed_out == -1 or self.rdkit_mol.GetNumConformers() == 0:
                print(f'2nd embed failed for {name} with SMILES: {smi}; will try standard embed')
                try:
                    embed_out = AllChem.EmbedMultipleConfs(self.rdkit_mol, numConfs=n_conformers)
                except:
                    embed_out = -1
                if embed_out == -1 or self.rdkit_mol.GetNumConformers() == 0:
                    print(f'3rd embed failed for {name} with SMILES: {smi}; wil try ETDG')
                    ps = AllChem.ETDG()
                    ps.randomSeed = 123
                    ps.useSmallRingTorsions=True
                    # ps.useMacrocycleTorsions=True
                    ps.ETversion=2
                    ps.useBasicKnowledge=True
                    try:
                        embed_out = AllChem.EmbedMultipleConfs(self.rdkit_mol, numConfs=n_conformers, params=ps)
                    except:
                        embed_out = -1
                    if embed_out == -1 or self.rdkit_mol.GetNumConformers() == 0:
                        raise Exception(f'4th embed failed for {name} with SMILES: {smi}')
        
        # Optimize structure with FF and get the lowest energy conformer
        energies = AllChem.MMFFOptimizeMoleculeConfs(self.rdkit_mol, maxIters=2000, nonBondedThresh=100.0)
        energies_list = [e[1] for e in energies]
        min_e_index = energies_list.index(min(energies_list))
        new_mol.AddConformer(self.rdkit_mol.GetConformer(min_e_index))
        self.rdkit_mol = new_mol

        # Make seperate directory for mol
        self.mol_calc_path = f'{self.SQMroot}/{name}'
        if not os.path.exists(self.mol_calc_path):
            os.mkdir(self.mol_calc_path)
        
        # Mol object to xyz file
        self.xyz_file_path = f'{self.mol_calc_path}/{name}.xyz'
        Chem.rdmolfiles.MolToXYZFile(self.rdkit_mol, self.xyz_file_path)


    def calc_CM5_charges(self, smi, name='pred_mol', optimize=False, save_output=False):
        """
        Run GFN1-xTB calculations to obtain CM5 atomic charges.
        :parameter optimize: if set to true, a GFN1-xTB geometry optimization is triggered.
        """

        # Generate xyz file from SMILES
        self.generate_3Dxyz(smi, name)

        # Get molecule properties
        chrg = Chem.GetFormalCharge(self.rdkit_mol)
        spin = 0 #spin hardcoded to zero

        # Run xTB calc
        if optimize:
            cmd = f'{XTBHOME}/bin/xtb --gfn 1 {self.xyz_file_path} --opt --lmo --chrg {chrg} --uhf {spin}' #TODO! add connectivity check!
        else:
            cmd = f'{XTBHOME}/bin/xtb --gfn 1 {self.xyz_file_path} --lmo --chrg {chrg} --uhf {spin}'
        
        proc = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, cwd=self.mol_calc_path)
        output = proc.communicate()[0]
        
        # Save calc output
        if save_output:
            with open(f'{self.mol_calc_path}/xtb.out', 'w') as f:
                f.write(output)

        # Get CM5 charges from output and append CM5 charges to RDKit mol object 
        self.cm5_list = []
        natoms = int(self.rdkit_mol.GetNumAtoms())
        for line_idx, line in enumerate(output.split('\n')):
            if 'Mulliken/CM5' in line:
                start = line_idx + 1
                endindex = start + natoms
                for i in range(start, endindex):
                    line = output.split('\n')[i]
                    cm5_atom = float(line.split()[2])
                    self.cm5_list.append(cm5_atom)
                break
        
        return self.cm5_list


    def create_descriptor_vector(self, atom_sites, n_shells=5, max_neighbors=4, use_cip_sort=True):
        """
        Create the GraphChargeShell descriptor for atoms in the list atom_sites.
        :parameter atom_sites example: [0,1]
        """

        NodeDescGen = NodeDescGenerator(molobj = self.rdkit_mol,
                                        property_list = self.cm5_list,
                                        n_shells = n_shells,
                                        max_neighbors = max_neighbors,
                                        use_cip_sort = use_cip_sort)

        # Create descriptor vector only for the provided atom sites     
        descriptor_vector = []
        mapper_vector = []
        for atom in atom_sites:
            atom_descriptor, atom_mapper = NodeDescGen.calculate_descriptor(atom)
            descriptor_vector.append(atom_descriptor)
            mapper_vector.append(atom_mapper)
        
        return descriptor_vector, mapper_vector
