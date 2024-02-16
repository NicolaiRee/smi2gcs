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

import logging
import numpy as np
import pandas as pd
from typing import Any
from rdkit import Chem

_logger = logging.getLogger(__name__)
#_logger.setLevel(logging.DEBUG)
_logger.setLevel(logging.WARNING)

class NodeDescGenerator():
    """Generate node descriptors based on 
       Cahn-Ingold-Prelog (CIP) or numerical sorting 
       of node features in a graph.
    """

    def __init__(
        self,
        molobj: object = None,
        property_list: list = None,
        n_shells: int = 5,
        max_neighbors: int = 4,
        use_cip_sort: bool = True,
        **kwargs: Any,
    ) -> None:
        
        super().__init__(**kwargs)

        self.molobj: object = molobj
        self.property_list: list = np.array(property_list)
        self.n_shells: int = n_shells
        self.max_neighbors: int = max_neighbors
        self.use_cip_sort: bool = use_cip_sort

        # The following can be replaced by a list of node types similiar to atomic numbers
        # and a custom adjacency matrix if RDKit molecule objects are not available.
        self.a_list = np.array([atom.GetAtomicNum() for atom in self.molobj.GetAtoms()]) # atomic numbers of atoms in molecule
        self.ac_matrix = Chem.rdmolops.GetAdjacencyMatrix(self.molobj, useBO=False) * self.a_list # adjacency matrix 

    
    def calculate_descriptor(self, atomIdx):

        descriptor_elements = [self.property_list[atomIdx]]
        mapper = [atomIdx]

        contained_atoms = set()  # store atom indices of all processed atoms
        contained_atoms.add(atomIdx)
        current_shell = [atomIdx]

        length = self.max_neighbors  # first atom has max_neighbors
        for shell in range(self.n_shells):
            
            last_shell = current_shell
            current_shell = []
            for ls_atom in last_shell:
                # Loop through atoms in shell and create next shell and sort block/subgraph
                if ls_atom != 'NaN':
                    # Get connected atoms
                    block = self._get_connected_atoms([[ls_atom]], contained_atoms)[0]
                    
                    # Sort the block according to the chosen sorting scheme    
                    if self.use_cip_sort:        
                        block = self._cip_sort_block(block, contained_atoms)
                    else:
                        block = self._sort_block(block)

                    # Update contained atoms
                    contained_atoms.update(block)
                else:
                    block = []

                # Fill block to length with 'NaN' if ls_atom do not have enough neighbors
                block = self._fill_block(block, length, item='NaN')

                # Update current shell
                current_shell.extend(block)

            descriptor_elements += [self.property_list[i] if i != 'NaN' else float(0.0) for i in current_shell]
            mapper += current_shell

            length = self.max_neighbors - 1  # one atoms is already contained, thus making length one less
        
        # Check if node descriptor has the correct length 
        assert len(descriptor_elements) == self._calculate_length(), f'{len(descriptor_elements)}  {self._calculate_length()}'
        
        return descriptor_elements, mapper
    

    def _calculate_length(self):
        result = 1 + self.max_neighbors
        previous_shell = self.max_neighbors

        for shell in range(1,self.n_shells):
            previous_shell = (self.max_neighbors - 1) * previous_shell
            result += previous_shell
        return result


    def _fill_block(self, block, length, item):
        """ Fills block to length with item """
        assert len(block) <= length, 'There is an atom with more than max neighbours!\n{}'.format(block)
        while len(block) < length:
            block.append(item)    
        return block


    def _get_connected_atoms(self, subgraphs, contained_atoms):
        """ Get idxs of connected atoms (subgraphs) """
        subgraphs = [np.nonzero(self.ac_matrix[i])[1] for i in subgraphs] # extract the next subgraphs
        subgraphs = [list(set(atoms) - set(contained_atoms)) for atoms in subgraphs] #exclude already included atoms
        return subgraphs
    

    def _sort_block(self, block):
        """ Sort the block according to the input properties in descending order.
        :param block: Block of atoms that need to be sorted (type: list)
        :return: sorted block (type: list)
        """
        block = sorted(block, key=lambda atom: float(self.property_list[atom]), reverse=True)
        return block


    def _cip_sort_block(self, block, contained_atoms): 
        """ Sort block according to CIP rules and input properties if CIP is unambiguous
    
            Sorting according to the CIP rules works as follows:
            1) Sort according to atomic number in descending order.
            2) If (1) is not unique, for each atom with the same priority (A*):
                a) Go to bound and yet not included atoms and sum up atomic numbers. Set the priority of A* according to summed atomic numbers.
                b) If (2a) did not give an unambiguous result expand the shell of each atom A* by one bond.
                c) repeat (2b) until a unique order is found.
            3) If no unique order is found in (2) and all bound atoms are included, then 
            sort atoms according to the input properties in descending order (this is an arbitrary choice).
        """

        _logger.debug(f'Entering CIP sorting routine: {block}')
        
        # Get priorities of atoms in shell
        priorities = self.a_list[block]

        # Create a dataframe with the priorities for the first bound atoms
        df_prior = pd.DataFrame(data=priorities, index=block, columns=[0])

        # Create a list of subgraphs
        subgraph = [[atom] for atom in block]

        while len(priorities) != len(set(priorities)):
            
            # NOTE ac_priorities could also be defined as np.sum(ac_matrix, axis=1) outside the while loop, 
            # but then subgraph atoms would be counted twice which makes it harder to debug, although the sorting would not change.
            not_included_atoms = [0 if i in contained_atoms else 1 for i in range(self.ac_matrix.shape[0])]
            ac_priorities = np.dot(self.ac_matrix, not_included_atoms)
            
            # Set priority of A* according to summed atomic numbers for bound and yet not included atoms
            # NOTE priorities += np.sum(ac_priorities[subgraphs], axis=1) would be an option,
            # if subgraphs had the same shape e.g. no subgraphs = [[13], [2], [10, 11, 27]]
            priorities = [int(np.sum(ac_priorities[s])) for s in subgraph] # converts the sum into integer as atomic numbers are integers

            # Break if subgraph has no further connections resulting in no changes to priorities
            if not any(priorities):
                break

            # Add subgraph priorities to df_prior
            df_prior[len(df_prior.columns)] = priorities
            
            # Remove all zero elements from priorities and concatenate previous priorities for the remainder
            priorities = df_prior[df_prior[df_prior.keys()[-1]] != 0].apply(lambda row: '-'.join(str(x) for x in row.tolist()), axis=1).tolist()
            _logger.debug(f'No unique order found for: {priorities}')

            # Update the list of included atoms
            contained_atoms = list(set(contained_atoms) | set(np.concatenate(subgraph)))
            
            # Get bound and yet not included atoms
            subgraph = self._get_connected_atoms(subgraph, contained_atoms)

        # Add input properties to df_prior
        df_prior['prop'] = np.array(self.property_list)[block]

        # Sort df_prior according to CIP rules and input properties if CIP is unambiguous
        sorted_df_prior = df_prior.sort_values(list(df_prior.keys()), ascending=False)
        _logger.debug(f'Sorted dataframe with the priorities: \n{sorted_df_prior}')
        
        sorted_block = sorted_df_prior.index.tolist()
        _logger.debug(f'Sorted block: {sorted_block}')
        
        return sorted_block