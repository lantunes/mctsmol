from mctsmol import *

if __name__ == '__main__':
    # allowed_angle_values = [0., 30., 60., 90., 120., 150., 180., -150., -120., -90., -60., -30.]
    allowed_angle_values = [0., 60., 120., 180., -120., -60.]
    num_simulations = 4000
    c = 3
    energy_function = mmff94_potential

    # mcts = MolecularMCTS3d(allowed_angle_values, energy_function, c=c)
    mcts = MolecularMCTS3e(allowed_angle_values, energy_function, energy_min=-5., energy_max=32., c=c)

    """
    - we can start with a randomly optimized structure, and just keep track of the lowest energy conformer
      found so far, and don't start a new search for each bond
    - we can see how many iterations we need to get a certain performance
      
      c=sqrt2, n=1000, 60 deg., molecule_96: -130.63581886662993
      
      c=sqrt2, n=20, 60 deg., molecule_2: -10.683416, -10.69823, -10.69823
      c=sqrt2, n=200, 60 deg., molecule_2: -10.69823, -10.69823, -10.69823
      c=sqrt2, n=1000, 60 deg., molecule_2: -10.69823, -10.69823
      c=sqrt2, n=3000, 30 deg., molecule_2: -10.69823
      
      c=5, n=1000, 60 deg., molecule_100: 1.8877847018712037
      c=0.5, n=1000, 60 deg., molecule_100: 1.2386237310353612
      c=sqrt2, n=1000, 60 deg., molecule_100: 0.6258685457052615, 0.4311348527216565
      c=sqrt2, n=2000, 60 deg., molecule_100: 0.4868096163827298
      c=sqrt2, n=3000, 60 deg., molecule_100: 0.18023328896400992
      c=sqrt2, n=2000, 30 deg., molecule_100: 0.11635177116116324, 2.0125215458601198
      c=sqrt2, n=3000, 30 deg., molecule_100: 0.615653774562932
      
      c=sqrt2, n=1000, 60 deg., molecule_71: 80.67395674652113, 80.55186999112318
      c=sqrt2, n=3000, 60 deg., molecule_71: 80.55186985680163
      
      c=sqrt2, n=1000, 60 deg., molecule_70: -43.00527190341086, -43.080350950899685
      
      c=sqrt2, n=1000, 60 deg., molecule_69: -3.4233533210070854, -4.640084095259487
      
      c=sqrt2, n=1000, 60 deg., molecule_98: 36.91644687543319
      
      for each molecule, c, n, deg.:
      - lowest energy found (10 trials)
      - mean lowest energy (10 trials)
      - std lowest energy (10 trials)
      - # of times lowest energy found (10 trials)
    """

    # state = mcts.init_state("FCCCCF")
    # state = mcts.init_state("c1(Cc2ccc(C#C[C@@H](N(C(N)=O)O)C)s2)ccc(F)cc1")  # molecule_2 (-10.69833, d=6)
    state = mcts.init_state("CCCCC[C@H](O)/C=C/[C@@H]1[C@@H](C/C=C\CCCC(O)=O)[C@@H]2C[C@H]1OO2")  # molecule_100 (-0.77875, d=14)
    # state = mcts.init_state("COC1=CC(N)=C(Cl)C=C1C(=O)N[C@H]1CCN(C[C@H]1OC)CCCOC1C=CC(F)=CC=1")  # molecule_71 (80.29538, d=11)
    # state = mcts.init_state("CC(C)/N=C(\\N)/N=C(\\N)/NOCCCOC1C=CC(=CC=1Cl)OC(F)(F)F")  # molecule_96 (-130.98873, d=13)
    # state = mcts.init_state("CCCCCNC(=N)N/N=C/C1=CNC2C=CC(=CC1=2)OC") # molecule_70 (-44.19546, d=9)
    # state = mcts.init_state("COC1C=C(CNC(=O)CCCC/C=C/C(C)C)C=CC=1O") # molecule_69 (-4.64008, d=11)
    # state = mcts.init_state("CC(C)C1N=C(C(C)C)C(COC)=C(C=1/C=C/[C@@H](O)C[C@@H](O)CC(O)=O)C1C=CC(F)=CC=1") # molecule_98 (36.40735, d=14)
    # state = mcts.init_state("CNCC[C@H](OC1C=CC(=CC=1)C(F)(F)F)C1C=CC=CC=1")  # molecule_57 (73.98046, d=7))

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
