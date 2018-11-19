from mctsmol import *

if __name__ == '__main__':
    num_angles = 10
    # allowed_angle_values = [60., -60., 120., -120., 0.]
    allowed_angle_values = [180., 60., -60.]
    num_simulations = 50
    # energy_function = simple_torsional_potential_minimized
    # energy_function = simple_torsional_potential
    energy_function = simple_torsional_potential_correlated

    mcts = TorsionalMCTS(num_angles, allowed_angle_values, energy_function, c=5)
    # mcts = TorsionalMCTS2(num_angles, allowed_angle_values, energy_function, c=5, initial_energy=5.)

    state = []

    while len(state) < num_angles:
        state = mcts.search(state, num_simulations)

    print("final state: %s (energy: %.5f)" %(state, energy_function(state)))
