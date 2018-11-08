from mctsmol import *

if __name__ == '__main__':
    trials = 100
    num_simulations = 300
    c = 100
    best_energies = []
    for i in range(0, trials):

        # allowed_angle_values = [0., 30., 60., 90., 120., 150., 180., -150., -120., -90., -60., -30.]
        allowed_angle_values = [0., 60., 120., 180., -120., -60.]

        energy_function = mmff94_potential

        mcts = MolecularMCTS3d(allowed_angle_values, energy_function, c=c)

        """
        molecule_57 (73.98046)
        497/279936 -> 73.98046 (0.18%)
          c=1, n=100:
            best energies mean: 74.06983
            best energies min: 73.98046
            global min times found: 11/30
          c=1, n=300:
            best energies mean: 73.99645
            best energies min: 73.98046
            global min times found: 21/30
            run #2:
            best energies mean: 74.00427
            best energies min: 73.98046
            global min times found: 71/100
          c=sqrt(2), n=300:
            best energies mean: 74.02309
            best energies min: 73.98046
            global min times found: 17/30
            run #2:
            best energies mean: 74.00673
            best energies min: 73.98046
            global min times found: 71/100
          c=2, n=300:
            best energies mean: 74.01357
            best energies min: 73.98046
            global min times found: 20/30
            run #2:
            best energies mean: 74.00604
            best energies min: 73.98046
            global min times found: 74/100
          c=2.5, n=300:
            best energies mean: 74.00977
            best energies min: 73.98046
            global min times found: 19/30
          c=3, n=300:
            best energies mean: 73.99949
            best energies min: 73.98046
            global min times found: 24/30
            run #2:
            best energies mean: 74.02271
            best energies min: 73.98046
            global min times found: 18/30
            run #3:
            best energies mean: 74.00489
            best energies min: 73.98046
            global min times found: 72/100
            run #4:
            best energies mean: 74.00992
            best energies min: 73.98046
            global min times found: 73/100
            run #5:
            best energies mean: 74.00729
            best energies min: 73.98046
            global min times found: 74/100
          c=3.5, n=300:
            best energies mean: 74.01395
            best energies min: 73.98046
            global min times found: 19/30
          c=4, n=300:
            best energies mean: 74.01357
            best energies min: 73.98046
            global min times found: 20/30
          c=5, n=100: 
            best energies mean: 74.08270
            best energies min: 73.98046
            global min times found: 10/30
          c=5, n=300:
            best energies mean: 74.02727
            best energies min: 73.98046
            global min times found: 17/30
          c=10, n=100:
            best energies mean: 74.09959
            best energies min: 73.98046
            global min times found: 12/30
          c=20, n=100:
            best energies mean: 74.08726
            best energies min: 73.98046
            global min times found: 11/30
          c=30, n=100:
            best energies mean: 74.07029
            best energies min: 73.98046
            global min times found: 11/30
          c=100, n=100:
            best energies mean: 74.12149
            best energies min: 73.98046
            global min times found: 14/30
          c=100, n=200:
            best energies mean: 74.01814
            best energies min: 73.98046
            global min times found: 19/30
          c=100, n=300:
            best energies mean: 73.98655
            best energies min: 73.98046
            global min times found: 25/30
            run #2:
            best energies mean: 74.00741
            best energies min: 73.98046
            global min times found: 73/100
          c=100, n=400:
            best energies mean: 74.00444
            best energies min: 73.98046
            global min times found: 22/30
          c=120, n=100:
            best energies mean: 74.04156
            best energies min: 73.98046
            global min times found: 14/30
          c=150, n=100:
            best energies mean: 74.08093
            best energies min: 73.98046
            global min times found: 15/30
          c=150, n=200:
            best energies mean: 74.03222
            best energies min: 73.98046
            global min times found: 15/30
          c=150, n=300:
            best energies mean: 74.00406
            best energies min: 73.98046
            global min times found: 23/30
          c=175, n=100:
            best energies mean: 74.07865
            best energies min: 73.98046
            global min times found: 15/30
          c=175, n=200:
            best energies mean: 74.01564
            best energies min: 73.98046
            global min times found: 21/30
          c=175, n=250:
            best energies mean: 74.01738
            best energies min: 73.98046
            global min times found: 21/30
          c=175, n=300:
            best energies mean: 74.01281
            best energies min: 73.98046
            global min times found: 22/30
          c=200, n=100:
            best energies mean: 74.09036
            best energies min: 73.98046
            global min times found: 9/30
          c=1000, n=100:
            best energies mean: 74.09837
            best energies min: 73.98046
            global min times found: 11/30
            
          --------  
          
          c=100, n=300, average value always 0:
            best energies mean: 74.01395
            best energies min: 73.98046
            global min times found: 19/30
            run #2:
            best energies mean: 74.01289
            best energies min: 73.98046
            global min times found: 69/100
          c=0, n=300:
            best energies mean: 74.17442
            best energies min: 73.98046
            global min times found: 8/30
          c=1, n=300, average value always 0:
            best energies mean: 74.00406
            best energies min: 73.98046
            global min times found: 23/30
            run #2:
            best energies mean: 74.00644
            best energies min: 73.98046
            global min times found: 72/100
          UCB1 formula with c=sqrt(2), n=300:
            best energies mean: 74.00375
            best energies min: 73.98046
            global min times found: 72/100
          UCB1 formula with c=sqrt(2), n=300, average value always 0:
            best energies mean: 74.01409
            best energies min: 73.98046
            global min times found: 71/100
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
