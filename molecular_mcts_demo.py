from mctsmol import *

if __name__ == '__main__':
    allowed_angle_values = [60., 120., 180., 240., 300.]
    num_simulations = 100
    energy_function = mmff94_potential

    mcts = MolecularMCTS(allowed_angle_values, energy_function, c=18)

    # state = mcts.init_state("C(C)(C)CCN")
    state = mcts.init_state("c1(Cc2ccc(C#C[C@@H](N(C(N)=O)O)C)s2)ccc(F)cc1")
    num_angles = mcts.get_num_angles()

    while len(state[1]) < num_angles:
        state = mcts.search(state, num_simulations)
        print("angles: %s; energy: %s" % (state[1], energy_function(state[0])))

    print("final state: %s (energy: %s)" %(state[1], energy_function(state[0])))

    with open('molecular_mcts.mol', 'w') as f:
        state[0].SetProp("_Name", "optimized_structure")
        f.write(Chem.MolToMolBlock(state[0]))
