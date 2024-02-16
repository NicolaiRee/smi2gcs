import copy
import logging
import numpy as np
from itertools import groupby

from DescriptorCreator.DescriptorElement import DescriptorElement


class GraphChargeShell(DescriptorElement):

    description = 'Charge shell based on sorted molecular graph'

    def __init__(self, **options):
        super(GraphChargeShell, self).__init__()
        self.n_shells = options.pop('n_shells')
        self.charge = options.pop('charge_type')
        self.cip = options.pop('use_cip_sort', False) # USE Cahn Ingold Prelog Type rules for sorting. If False, then sort purely based on charge
        self._current_mol = None

    def calculate_elements(self, atom):
        charge = self.charge
        self._current_mol = atom.GetOwningMol()
        last_shell = []
        current_shell = []
        block = []  # Corresponds to a block of atoms which then have to be sorted, could be sourced to a own class

        # set the atomic charge as the first element
        descriptor_elements = [float(atom.GetProp(charge))]
        mapper = [atom.GetIdx()]
        contained_atoms = set()  # Store atom indices of all processed atoms
        contained_atoms.add(atom.GetIdx())
        last_shell.append(atom)

        length = 4  # first atom has 4 neighbours
        for shell in range(self.n_shells):
            for ls_atom in last_shell:
                # loop through shell and create next shell
                # create and sort subshells / blocks
                for neighbour_atom in ls_atom.GetNeighbors():
                    if neighbour_atom.GetIdx() not in contained_atoms:
                        block.append(neighbour_atom)

                # Sort the block according to the chosen sorting scheme
                if len(block):
                    if self.cip:
                        block = self._sort_block_cip(block, contained_atoms)
                    else:
                        block = self._sort_block(block)
                
                # Fill the block with dummy atoms, such that all blocks have the same length
                self._fill_block(block, length)

                current_shell = current_shell + block
                block = []
                length = 3  # from the second shell on atoms only have 3 members because one is already contained

                indices_block = [atm.GetIdx() for atm in current_shell]
                contained_atoms.update(indices_block) 

            indices = [atm.GetIdx() for atm in current_shell]
            descriptor_elements += [float(atm.GetProp(charge)) for atm in current_shell]
            mapper += indices
            last_shell = current_shell
            current_shell = []
        assert len(descriptor_elements) == self._calculate_length(self.n_shells), '{}  {}'.format(len(
            descriptor_elements),  self._calculate_length(self.n_shells))  # Debugging assertion
        
        return descriptor_elements, mapper


    @classmethod
    def _calculate_length(cls, n):
        result = 5
        previous_shell = 4
        m = 1
        while m < n:
            result += 3*previous_shell
            previous_shell = 3*previous_shell
            m += 1
        return result


    def _fill_block(self, block, length):
        """ Fills block with dummy atoms to length and than further sends it to the sorting routine
        :param block: block containing the atoms to sort
        :param length: the length to which the block should be filled with dummy atoms
        :return:
        """

        class DummyA(object):
            def __init__(self):
                self.charge = 0.0
            
            @classmethod
            def GetNeighbors(cls):
                return []

            @classmethod
            def GetIdx(cls):
                return -1

            def GetProp(self, charge):
                return self.charge

        assert len(block) <= length, 'There is a atom with more than 4 neighbours!\n{}'.format(block)
        while len(block) < length:
            block.append(DummyA())


    def _sort_block(self, block):
        """ Sort the block according to charges in descending order.
        :param block: Block of Atoms that need to be sorted (type: list)
        :type block: list
        :return: sorted block (type: list)
        """
        block = sorted(block, key=lambda atom: float(atom.GetProp(self.charge)), reverse=True)
        return block
    

    def _expand_shell(self, atom_list, contained_atoms):
        subgraph_temp = []
        priorities = 0

        for atom in atom_list:
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in contained_atoms:
                    subgraph_temp.append(nbr)
                    priorities += nbr.GetAtomicNum()

        return subgraph_temp, priorities


    def _sort_block_cip(self, block, contained_atoms):
        """ Sort block according to CIP rules and charges if CIP is unambiguous

        Sorting according to the CIP rules works as follows:
        1) Take all bound substituents and sort according to atomic number in descending order.
        2) If 1 is not unique, for each atom with same priority (A*):
            a) Go to bound and yet not included atoms and sum up atomic numbers. Set priority of A* according to summed atomic numbers.
            b) If 2 a) did not give unambiguous result expand shell of each atom A* by one bond.
            c) repeat 2b) until unique order is found.
            d) If no unique order is found and all bound atoms are included (in set of already summed in atoms)
                sort atoms according to charge in descending order (this is an arbitrary choice).

        :param block: Block of Atoms that need to be sorted
        :type block: list
        :return: CIP sorted block (type: list)
        """
        
        logging.debug('Entering CIP determination routine')
        logging.debug(f'Length of block: {len(block)}')

        contained_atoms_in_sort = set(contained_atoms)
        contained_atoms_in_sort.update(atom.GetIdx() for atom in block)
        logging.debug(f'Contained atoms in sort: {contained_atoms_in_sort}')

        subgraphs = [[atm] for atm in block]
        priorities = [atm.GetAtomicNum() for atm in block]

        while len(priorities) != len(set(priorities)):
            old_priorities = copy.deepcopy(priorities)

            for num, atom_list in enumerate(subgraphs):
                subgraph_temp, priority_temp = self._expand_shell(atom_list, contained_atoms_in_sort)
                subgraphs[num] = subgraph_temp
                priorities[num] += priority_temp

            contained_atoms_in_sort.update(atom.GetIdx() for atom in np.concatenate(subgraphs))

            if priorities == old_priorities:
                logging.debug(f'No unique order found: {priorities}, {list(set(priorities))}')
                break

        # Sort block according the priorities
        priorities, block = zip(*sorted(zip(priorities, block), key=lambda x:x[0], reverse=True))
        priorities, block = list(priorities), list(block)
        logging.debug(f'Priorities after sort: {priorities}')
        
        # Now check if all the priorities are different and if not 
        # sort those with the same priorities according to the charge in descending order
        if len(priorities) != len(set(priorities)):
            logging.debug(f'Charges before sorting: {[atom.GetProp(self.charge) for atom in block]}')

            # Group atoms with the same priority p=[(8,atom1),(8,atom2),(1,atom3)] --> [[(8,atom1),(8,atom2)], [(1,atom3)]]
            sets_to_sort = [list(j) for i, j in groupby(zip(priorities,block), key=lambda x:x[0])]
            logging.debug(f'Grouped atoms (idx): {[list(np.array(list(j))[:,1]) for i, j in groupby(zip(priorities,[atom.GetIdx() for atom in block]), key=lambda x:x[0])]}')

            # Sort grouped atoms according to the charge in descending order 
            sorted_block = []
            for index_set in sets_to_sort:
                sub_block = [x[1] for x in index_set]
                sorted_block.extend(sorted(sub_block, key=lambda atom:float(atom.GetProp(self.charge)), reverse=True))
            
            assert len(sorted_block) == len(block)
            block = sorted_block
            logging.debug(f'Charges after sorting: {[atom.GetProp(self.charge) for atom in block]}')

        return block
