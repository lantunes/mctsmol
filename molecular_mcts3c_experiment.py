from mctsmol import *

if __name__ == '__main__':
    trials = 100
    num_simulations = 300
    c = 50
    best_energies = []
    for i in range(0, trials):

        # allowed_angle_values = [0., 30., 60., 90., 120., 150., 180., -150., -120., -90., -60., -30.]
        allowed_angle_values = [0., 60., 120., 180., -120., -60.]

        energy_function = mmff94_potential

        mcts = MolecularMCTS3c(allowed_angle_values, energy_function, c=c)

        """
        molecule_57 (73.98046)
        497/279936 -> 73.98046 (0.18%)
          c=50, n=300, UCB1:
            best energies mean: 73.99507
            best energies min: 73.98046
            global min times found: 82/100
            run #2:
            best energies mean: 73.99964
            best energies min: 73.98046
            global min times found: 75/100
            run #3:
            best energies mean: 74.00615
            best energies min: 73.98046
            global min times found: 73/100
            run #4:
            best energies mean: 74.01312
            best energies min: 73.98046
            global min times found: 70/100
          c=5, n=300, UCB1:
            best energies mean: 74.01026
            best energies min: 73.98046
            global min times found: 70/100
          c=100, n=300, UCB1:
            best energies mean: 73.99720
            best energies min: 73.98046
            global min times found: 78/100
          c=35, n=300, UCB1:
            best energies mean: 74.00432
            best energies min: 73.98046
            global min times found: 78/100
            run #2:
            best energies mean: 74.01044
            best energies min: 73.98046
            global min times found: 70/100
          c=150, n=300, UCB1:
            best energies mean: 74.01181
            best energies min: 73.98046
            global min times found: 69/100
            
          c=5, n=300, PUCT:
            best energies mean: 74.04785
            best energies min: 73.98046
            global min times found: 55/100
          c=35, n=300, PUCT:
            best energies mean: 74.01060
            best energies min: 73.98046
            global min times found: 67/100
          c=75, n=300, PUCT:
            best energies mean: 74.01283
            best energies min: 73.98046
            global min times found: 71/100
          c=100, n=300, PUCT:
            best energies mean: 74.01437
            best energies min: 73.98046
            global min times found: 67/100
        """

        # state = mcts.init_state("FCCCCF")
        # state = mcts.init_state("c1(Cc2ccc(C#C[C@@H](N(C(N)=O)O)C)s2)ccc(F)cc1")  # molecule_2 (-10.69833, d=6)
        state = mcts.init_state("CNCC[C@H](OC1C=CC(=CC=1)C(F)(F)F)C1C=CC=CC=1")  # molecule_57 (73.98046, d=7))

        mcts.search(state, num_simulations)

        best_conformer = mcts.get_best_conformer()
        best_conformer_energy = mcts.get_best_conformer_energy()
        dihedrals = mcts.get_dihedrals(best_conformer)
        print("best conformer energy: %s" % best_conformer_energy)
        print("best conformer dihedrals: %s" % dihedrals)
        print("best rollout state: %s" % mcts.get_best_rollout_state())

        with open('molecular_mcts3.mol', 'w') as f:
            best_conformer.SetProp("_Name", "optimized_structure")
            f.write(Chem.MolToMolBlock(best_conformer))

        best_energies.append(best_conformer_energy)

    print("--------------------")
    print("c=%s, n=%s:" % (c, num_simulations))
    print("best energies: %s" % ["%.5f" % x for x in best_energies])
    print("best energies mean: %.5f" % np.mean(best_energies))
    print("best energies min: %.5f" % np.min(best_energies))
    unique, counts = np.unique(["%.5f" % x for x in best_energies], return_counts=True)
    d = dict(zip(unique, counts))
    print("global min times found: %s/%s" % (d['73.98046'], trials))
