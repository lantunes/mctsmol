from mctsmol import *

if __name__ == '__main__':
    num_angles = 10
    allowed_angle_values = [180., 60., -60.]
    num_simulations = 50

    mcts = TorsionalMCTS(num_angles, allowed_angle_values, simple_torsional_potential, c=5)

    state = []

    while len(state) < num_angles:
        state = mcts.search(state, num_simulations)

    print("final state: %s (energy: %s)" %(state, simple_torsional_potential(state)))
