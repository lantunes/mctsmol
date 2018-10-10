import numpy as np


def simple_torsional_potential(angles):
    """
    :param angles: a list of angles in degrees 
    :return: the total torsional potential
    """
    sum = 0.0
    for angle in angles:
        sum += 1 + np.cos(np.deg2rad(angle)) + np.cos(np.deg2rad(3 * angle))
    return sum

# print(simple_torsional_potential([180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0]))