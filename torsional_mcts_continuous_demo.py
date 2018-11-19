from mctsmol import *

if __name__ == '__main__':
    allowed_angle_values = [60., 180., -60.]
    # allowed_angle_values = [0., 60., 120., 180., -60., -120.]
    # allowed_angle_values = [0., 60., 120., -60., -120.]
    num_angles = 10
    num_simulations = 150

    energy_function = simple_torsional_potential
    # energy_function = simple_torsional_potential_correlated
    # energy_function = simple_torsional_potential_minimized

    # mcts = TorsionalMCTSContinuous(num_angles, allowed_angle_values, energy_function, c=sqrt(2), initial_energy=5.)
    mcts = TorsionalMCTSContinuous2(num_angles, allowed_angle_values, energy_function, c=5)
    # mcts = TorsionalMCTSContinuous3(num_angles, allowed_angle_values, energy_function, c=5)

    state = []

    mcts.search(state, num_simulations)

    best_conformer = mcts.get_best_conformer()
    best_conformer_energy = mcts.get_best_conformer_energy()
    print("best conformer energy: %.5f" % best_conformer_energy)
    print("best rollout state: %s" % best_conformer)
