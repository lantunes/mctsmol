import numpy as np
from math import cos
from math import sin


def simple_torsional_potential_minimized(angles):
    """
    :param angles: a list of angles in degrees 
    :return: a tuple of the total minimized torsional potential and the dihedral angles after minimization 
    """
    rads = [np.deg2rad(x) for x in angles]
    res = gradient_descent(rads, partial_derivative, partial_second_derivative)
    return energy(res), np.rad2deg(res)


def energy(angles_in_rads):
    sum = 0.0
    for angle in angles_in_rads:
        sum += 1 + cos(angle) + cos(3 * angle)
    return sum


def partial_derivative(rads):
    return -3*sin(3*rads) - sin(rads)


def partial_second_derivative(rads):
    return -9*cos(3*rads) - cos(rads)


def gradient_descent(x, dfdx, d2dx2, learning_rate=0.001, convergence_threshold=1e-5, debug=False):
    x_new = list(x)
    grad = []
    grad2 = []
    for x_i in x:
        grad.append(dfdx(x_i))
        grad2.append(d2dx2(x_i))

    def is_converged(partial_derivatives, partial_second_derivates):
        for i,g in enumerate(partial_derivatives):
            g2 = partial_second_derivates[i]
            # if the second derivative is negative and the gradient is 0., then we are at a maxima
            if np.abs(g) > convergence_threshold or (g2 < 0 and g == 0.) :
                return False
        return True

    while not is_converged(grad, grad2):

        grad = []
        grad2 = []
        for x_i in x_new:
            g = dfdx(x_i)
            grad.append(-0.01 if g == 0. else g)
            grad2.append(d2dx2(x_i))

        for i,_ in enumerate(x_new):
            x_new[i] = x_new[i] - learning_rate*grad[i]

        if debug: print("grad: %s; x: %s" % (grad, x_new))

    return x_new


if __name__ == '__main__':
    # d/dx of energy: −3sin(3x)−sin(x)

    angles = [60., 120., -60., -120., 180., 60., 0., -60., 120., 120.]
    rads = [np.deg2rad(x) for x in angles]
    print("energy: %s" % energy(rads))
    print("rads: %s" % rads)

    res = gradient_descent(rads, partial_derivative, partial_second_derivative, debug=True)

    print("degrees: %s " % np.rad2deg(res))
    print("final energy: %s" % energy(res))
