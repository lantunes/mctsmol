import random
from math import *
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
from rdkit.Chem import PDBWriter
from rdkit.Chem import MolFromPDBFile
import copy

"""
Continuous MCTS. Reward squashed between -1 and 1.
"""
class MolecularMCTS3d:
    def __init__(self, allowed_angle_values, energy_function, c=sqrt(2)):
        self._allowed_angle_values = allowed_angle_values
        self._energy_function = energy_function
        self._c = c
        self._original_mol = None
        self._global_minimum_mol = None
        self._global_minimum_energy = None
        self._global_minimum_rollout_state = None

    def init_state(self, smiles_string):
        mol = Chem.MolFromSmiles(smiles_string)
        mol = Chem.AddHs(mol,explicitOnly=False)
        AllChem.EmbedMolecule(mol)

        mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp)
        opt_fail = ff.Minimize(maxIts=1000,forceTol=0.0001,energyTol=1e-06)
        energy = ff.CalcEnergy()

        self._rotatable_bonds = self._get_rotatable_bonds(mol)
        self._num_angles = len(self._rotatable_bonds)

        print("initial MMFF94 energy: %s" % energy)
        print("initial geometry: %s" % self.get_dihedrals(mol))

        self._original_smiles = self._get_smiles(mol)
        self._original_mol = mol
        return []

    def _get_rotatable_bonds(self, molecule):
        raw_rot_bonds =  molecule.GetSubstructMatches(Chem.MolFromSmarts("[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]"))
        raw_rot_bonds += molecule.GetSubstructMatches(Chem.MolFromSmarts("[*]~[*]-[O,S]-[#1]"))
        raw_rot_bonds += molecule.GetSubstructMatches(Chem.MolFromSmarts("[*]~[*]-[NX3;H2]-[#1]"))
        bonds = []
        rot_bonds = []
        for k,i,j,l in raw_rot_bonds:
            if (i,j) not in bonds and (j,i) not in bonds: # makes sure that dihedrals are unique
                bonds.append((i,j))
                rot_bonds.append((k,i,j,l))

        return rot_bonds

    def get_num_angles(self):
        if self._num_angles is None:
            raise Exception("init_state must be called first")
        return self._num_angles

    def search(self, state, num_simulations):
        """
        :param state: list of selected dihedral angles 
        :param num_simulations: the number of simulations to perform
        """
        root_node = _Node(state, self._rotatable_bonds, self._allowed_angle_values, self._c)

        # Perform simulations
        for i in range(num_simulations):
            node = root_node

            # Select
            while not node.has_untried_moves() and node.has_children():
                node = node.select_child()

            # Expand
            if node.has_untried_moves():
                move_state = node.select_untried_move()
                node = node.add_child(move_state, self._rotatable_bonds, self._allowed_angle_values, self._c)

            # Rollout
            rollout_state = node.state
            while len(rollout_state) < self._num_angles:
                rollout_state = self._select_next_move_randomly(rollout_state)

            # Backpropagate
            #   backpropagate from the expanded node and work back to the root node
            mol = self.get_mol_with_dihedrals(rollout_state)
            energy = self._energy_function(mol)
            reward = self._get_reward(energy)
            is_valid = self._is_valid_structure(mol)
            print("energy: %s; reward: %s; is valid: %s" % (energy, reward, is_valid))
            if (self._global_minimum_energy is None or energy < self._global_minimum_energy) and is_valid:
                self._global_minimum_energy = energy
                self._global_minimum_mol = mol
                self._global_minimum_rollout_state = rollout_state
            while node is not None:
                node.visits += 1
                node.rewards.append(reward if is_valid else -1.0)
                node = node.parent

    def _get_reward(self, energy):
        return -(energy / (1 + np.abs(energy)))

    def _is_valid_structure(self, mol):
        # checks if the stereochemistry of the molecule is conserved
        try:
            # somekind of bug here. Making the smiles from the mol object sometimes result in a wrong
            # stereochemistry, but writing to a pdb file and the reading it produces the correct result
            # this is a hack that should be fixed
            temp_pdb_name = 'isvalid_temp.pdb'
            self._write_pdb(mol, temp_pdb_name)
            temp_mol = MolFromPDBFile(temp_pdb_name)
            test_smiles = self._get_smiles(temp_mol)
            if self._original_smiles != test_smiles:
                return False
        except:
            # this is need to prevent rdkit from crashing if it cant understand the molecule
            print("WARNING: can't understand molecule")
            return False
        return True

    def _write_pdb(self, molecule, filename):
        w = PDBWriter('./'+filename)
        w.write(molecule)

    def _get_smiles(self, molecule):
        molecule_copy = copy.deepcopy(molecule)
        test_molecule = Chem.Mol(molecule_copy)
        test_molecule = Chem.RemoveHs(test_molecule)
        Chem.SanitizeMol(test_molecule)
        Chem.DetectBondStereochemistry(test_molecule,-1)
        Chem.AssignStereochemistry(test_molecule, flagPossibleStereoCenters=True, force=True)
        Chem.AssignAtomChiralTagsFromStructure(test_molecule,-1)
        smiles = Chem.MolToSmiles(test_molecule, isomericSmiles=True)
        return smiles

    def _select_next_move_randomly(self, state):
        return list(state) + [np.random.choice(self._allowed_angle_values)]

    def get_mol_with_dihedrals(self, state):
        molecule_copy = copy.deepcopy(self._original_mol)
        mol = Chem.Mol(molecule_copy)  # create a copy of the Mol
        conf = mol.GetConformer()
        for i, bond in enumerate(self._rotatable_bonds):
            rdMolTransforms.SetDihedralDeg(conf, bond[0], bond[1], bond[2], bond[3], state[i])
        conformer_id = mol.AddConformer(conf, assignId=True)
        return Chem.Mol(mol, False, conformer_id)

    def get_dihedrals(self, mol):
        conf = mol.GetConformer()
        dihedrals = []
        for bond in self._rotatable_bonds:
            dihedrals.append(rdMolTransforms.GetDihedralDeg(conf, bond[0], bond[1], bond[2], bond[3]))
        return dihedrals

    def get_best_conformer(self):
        return self._global_minimum_mol

    def get_best_conformer_energy(self):
        return self._global_minimum_energy

    def get_best_rollout_state(self):
        return self._global_minimum_rollout_state


class _Node:
    def __init__(self, state, rotatable_bonds, allowed_angle_values, c, parent=None):
        self.state = state
        self._rotatable_bonds = rotatable_bonds
        self._c = c
        self._allowed_angle_values = allowed_angle_values
        self.rewards = []
        self.visits = 0.0
        self.parent = parent
        self.children = []
        self.untried_moves = self._get_child_states()

    def _get_child_states(self):
        child_states = []
        if len(self.state) < len(self._rotatable_bonds):
            for allowed_angle_value in self._allowed_angle_values:
                child_states.append(list(self.state) + [allowed_angle_value])
        return child_states

    def _average_value(self):
        return np.sum(self.rewards) / self.visits

    def has_untried_moves(self):
        return self.untried_moves != []

    def select_untried_move(self):
        return random.choice(self.untried_moves)

    def add_child(self, child_state, rotatable_bonds, allowed_angle_values, c):
        child = _Node(child_state, rotatable_bonds, allowed_angle_values, c, parent=self)
        self.children.append(child)
        self.untried_moves.remove(child_state)
        return child

    def has_children(self):
        return self.children != []

    def select_child(self):
        highest_ucb1 = None
        selected_child_node = None
        for child_node in self.children:
            ucb1 = child_node.ucb1()
            if highest_ucb1 is None or highest_ucb1 < ucb1:
                highest_ucb1 = ucb1
                selected_child_node = child_node
        return selected_child_node

    def ucb1(self):
        if self.visits == 0:
            return math.inf

        prob = 1.0 / len(self._allowed_angle_values)  # equal probability for visiting each child node
        return self._average_value() + self._c*prob*(sqrt(self.parent.visits)/(1 + self.visits))

        # return self._average_value() + sqrt(2)*sqrt(log(self.parent.visits)/self.visits)
