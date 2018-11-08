from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms

from mctsmol import *


def init_state(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    mol = Chem.AddHs(mol,explicitOnly=False)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    return mol

def get_mol_with_dihedrals(original_mol, rotatable_bonds, state):
    mol = Chem.Mol(original_mol)  # create a copy of the Mol
    conf = mol.GetConformer()
    for i, bond in enumerate(rotatable_bonds):
        rdMolTransforms.SetDihedralDeg(conf, bond[0], bond[1], bond[2], bond[3], state[i])
    conformer_id = mol.AddConformer(conf, assignId=True)
    return Chem.Mol(mol, False, conformer_id)


def get_rotatable_bonds(molecule):
    raw_rot_bonds =  molecule.GetSubstructMatches(Chem.MolFromSmarts("[!#1]~[!$(*#*)&!D1]-!@[!$(*#*)&!D1]~[!#1]"))
    raw_rot_bonds += molecule.GetSubstructMatches(Chem.MolFromSmarts("[*]~[*]-[O,S]-[#1]"))
    raw_rot_bonds += molecule.GetSubstructMatches(Chem.MolFromSmarts("[*]~[*]-[NX3;H2]-[#1]"))
    bonds = []
    rot_bonds = []
    for k, i, j, l in raw_rot_bonds:
        if (i, j) not in bonds:
            bonds.append((i,j))
            rot_bonds.append((k, i, j, l))
    return rot_bonds


def get_dihedrals(mol, rotatable_bonds):
    conf = mol.GetConformer()
    dihedrals = []
    for bond in rotatable_bonds:
        dihedrals.append(rdMolTransforms.GetDihedralDeg(conf, bond[0], bond[1], bond[2], bond[3]))
    return dihedrals


if __name__ == '__main__':

    energy_function = mmff94_potential

    mol = init_state("CCCCC[C@H](O)/C=C/[C@@H]1[C@@H](C/C=C\CCCC(O)=O)[C@@H]2C[C@H]1OO2")
    dihedrals = [60., 60., 60., 60., 60., 60., 60., 60., 60., 60., 60., 60., 60., 60., 60.]

    bonds = get_rotatable_bonds(mol)

    print(bonds)

    mol = get_mol_with_dihedrals(mol, bonds, dihedrals)

    energy = energy_function(mol)

    print(energy)
    print(get_dihedrals(mol, bonds))



