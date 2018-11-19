from mctsmol import *

if __name__ == '__main__':
    allowed_angle_values = [60., 120., 180., 240., 300.]
    num_simulations = 100
    c = 5
    # energy_min = 50.
    # energy_max = 100.
    energy_min = -5.
    energy_max = 50.
    energy_function = mmff94_potential

    mcts = MolecularMCTS1b(allowed_angle_values, energy_function, energy_min=energy_min, energy_max=energy_max, c=c)

    # state = mcts.init_state("c1(Cc2ccc(C#C[C@@H](N(C(N)=O)O)C)s2)ccc(F)cc1")  # molecule_2
    state = mcts.init_state("CCCCC[C@H](O)/C=C/[C@@H]1[C@@H](C/C=C\CCCC(O)=O)[C@@H]2C[C@H]1OO2")  # molecule_100
    # state = mcts.init_state("COC1=CC(N)=C(Cl)C=C1C(=O)N[C@H]1CCN(C[C@H]1OC)CCCOC1C=CC(F)=CC=1")  # molecule_71
    # state = mcts.init_state("CC(C)/N=C(\\N)/N=C(\\N)/NOCCCOC1C=CC(=CC=1Cl)OC(F)(F)F")  # molecule_96
    # state = mcts.init_state("CNCC[C@H](OC1C=CC(=CC=1)C(F)(F)F)C1C=CC=CC=1")  # molecule_57 (73.98046, d=7))

    num_angles = mcts.get_num_angles()

    while len(state[1]) < num_angles:
        state = mcts.search(state, num_simulations)
        print("dihedrals: %s" % state[1])

    mol = mcts.get_mol_with_dihedrals(state[1])
    final_energy = energy_function(mol)
    final_dihedrals = mcts.get_dihedrals(mol)
    print("final state: %s (energy: %s)" %(final_dihedrals, final_energy))

    with open('molecular_mcts.mol', 'w') as f:
        state[0].SetProp("_Name", "optimized_structure")
        f.write(Chem.MolToMolBlock(state[0]))
