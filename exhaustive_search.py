import itertools

from mctsmol import *


def init_state(smiles_string):
    mol = Chem.MolFromSmiles(smiles_string)
    mol = Chem.AddHs(mol,explicitOnly=False)
    AllChem.EmbedMolecule(mol)

    write_pdb(mol, 'exhaustive_preoptimised.pdb')

    AllChem.MMFFOptimizeMolecule(mol)
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp)

    Chem.SanitizeMol(mol)
    Chem.DetectBondStereochemistry(mol,-1)
    Chem.AssignStereochemistry(mol, flagPossibleStereoCenters=True, force=True)
    Chem.AssignAtomChiralTagsFromStructure(mol,-1)

    opt_fail = ff.Minimize(maxIts=1000,forceTol=0.0001,energyTol=1e-06)
    energy = ff.CalcEnergy()

    write_pdb(mol, 'exhaustive_optimised.pdb')

    return mol, energy


def get_mol_with_dihedrals(original_mol, rotatable_bonds, state):
    molecule_copy = copy.deepcopy(original_mol)
    mol = Chem.Mol(molecule_copy)  # create a copy of the Mol
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
    for k,i,j,l in raw_rot_bonds:
        if (i,j) not in bonds and (j,i) not in bonds: # makes sure that dihedrals are unique
            bonds.append((i,j))
            rot_bonds.append((k,i,j,l))

    return rot_bonds


def get_dihedrals(mol, rotatable_bonds):
    conf = mol.GetConformer()
    dihedrals = []
    for bond in rotatable_bonds:
        dihedrals.append(rdMolTransforms.GetDihedralDeg(conf, bond[0], bond[1], bond[2], bond[3]))
    return dihedrals


def is_valid_structure(mol, original_smiles):
    # checks if the stereochemistry of the molecule is conserved
    try:
        # somekind of bug here. Making the smiles from the mol object sometimes result in a wrong
        # stereochemistry, but writing to a pdb file and the reading it produces the correct result
        # this is a hack that should be fixed
        temp_pdb_name = 'isvalid_temp.pdb'
        write_pdb(mol, temp_pdb_name)
        temp_mol = MolFromPDBFile(temp_pdb_name)
        test_smiles = get_smiles(temp_mol)
        if original_smiles != test_smiles:
            return False
    except:
        # this is need to prevent rdkit from crashing if it cant understand the molecule
        print("WARNING: can't understand molecule")
        return False
    return True


def write_pdb(molecule, filename):
    w = PDBWriter('./'+filename)
    w.write(molecule)


def get_smiles(molecule):
    molecule_copy = copy.deepcopy(molecule)
    test_molecule = Chem.Mol(molecule_copy)
    test_molecule = Chem.RemoveHs(test_molecule)
    Chem.SanitizeMol(test_molecule)
    Chem.DetectBondStereochemistry(test_molecule,-1)
    Chem.AssignStereochemistry(test_molecule, flagPossibleStereoCenters=True, force=True)
    Chem.AssignAtomChiralTagsFromStructure(test_molecule,-1)
    smiles = Chem.MolToSmiles(test_molecule, isomericSmiles=True)
    return smiles

if __name__ == '__main__':
    energy_function = mmff94_potential
    allowed_angle_values = [0., 60., 120., 180., -120., -60.]

    # smiles = "FCCCCF"  #(d=3)
    # smiles = "c1(Cc2ccc(C#C[C@@H](N(C(N)=O)O)C)s2)ccc(F)cc1"  # molecule_2 (d=6)
    smiles = "CNCC[C@H](OC1C=CC(=CC=1)C(F)(F)F)C1C=CC=CC=1"  # molecule_57 (d=7)

    mol, initial_energy = init_state(smiles)
    original_smiles = get_smiles(mol)

    rotatable_bonds = get_rotatable_bonds(mol)
    number_of_bonds = len(rotatable_bonds)

    print("initial MMFF94 energy: %s" % initial_energy)
    print("initial geometry: %s" % get_dihedrals(mol, rotatable_bonds))

    dihedral_permutations = list(itertools.product(allowed_angle_values, repeat=number_of_bonds))

    print("number of dihedral permutations: %s" % len(dihedral_permutations))

    energies = {}
    for dihedral_setting in dihedral_permutations:
        mol_with_dihedrals = get_mol_with_dihedrals(mol, rotatable_bonds, dihedral_setting)
        energy = "%.5f" % energy_function(mol_with_dihedrals)  # optimize the geometry
        is_valid = is_valid_structure(mol_with_dihedrals, original_smiles)
        actual_dihedrals = ["%.5f" % x for x in get_dihedrals(mol_with_dihedrals, rotatable_bonds)]
        print("%s; %s; %s; %s" % (dihedral_setting, energy, actual_dihedrals, is_valid))
        if is_valid:
            energies[energy] = energies.get(energy, 0) + 1

    print("number of unique minima: %s" % len(energies))
    print("energies by lowest: %s" % list(sorted(energies.items(), key=lambda kv: float(kv[0]))))
    print("energies by most frequent: %s" % list(reversed(sorted(energies.items(), key=lambda kv: kv[1]))))
