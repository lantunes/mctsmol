from mctsmol import *

if __name__ == '__main__':
    allowed_angle_values = [0., 60., 120., 180., -120., -60.]
    # allowed_angle_values = [0., 30., 60., 90., 120., 150., 180., -150., -120., -90., -60., -30.]
    energy_function = mmff94_potential
    num_simulations = 100
    c = 100

    mcts = MolecularMCTS2(allowed_angle_values, energy_function= energy_function, c=c)

    # state = mcts.init_state("c1(Cc2ccc(C#C[C@@H](N(C(N)=O)O)C)s2)ccc(F)cc1")  # molecule_2 (-10.69833)
    # state = mcts.init_state("CCCCC[C@H](O)/C=C/[C@@H]1[C@@H](C/C=C\CCCC(O)=O)[C@@H]2C[C@H]1OO2")  # molecule_100 (-0.77875)
    # state = mcts.init_state("COC1=CC(N)=C(Cl)C=C1C(=O)N[C@H]1CCN(C[C@H]1OC)CCCOC1C=CC(F)=CC=1")  # molecule_71 (80.29538)
    state = mcts.init_state("CC(C)/N=C(\\N)/N=C(\\N)/NOCCCOC1C=CC(=CC=1Cl)OC(F)(F)F")  # molecule_96 (-130.98873)

    num_angles = mcts.get_num_angles()

    while len(state) < num_angles:
        state = mcts.search(state, num_simulations)
        print("%s" % state)

    mol = mcts.get_mol_with_dihedrals(state)
    final_energy = energy_function(mol)
    final_dihedrals = mcts.get_dihedrals(mol)
    print("final state: %s (energy: %s)" %(final_dihedrals, final_energy))

    with open('molecular_mcts2.mol', 'w') as f:
        mol.SetProp("_Name", "optimized_structure")
        f.write(Chem.MolToMolBlock(mol))
