#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt


class Model:
    """Stochastic model constants."""
    mu = 3
    sigma = 1


def dW(dt):
    """Random sample normal distribution."""
    return np.random.normal(loc=0.0, scale=np.sqrt(dt))


def run_simulation():
    """ Return the result of one full simulation."""
    # One second and thousand grid points
    T_INIT = 0
    T_END = 1
    N = 1000  # Compute 1000 grid points
    DT = float(T_END - T_INIT) / N
    TS = np.arange(T_INIT, T_END + DT, DT)

    Y_INIT = 1

    # Vectors to fill
    ys = np.zeros(N + 1)
    ys[0] = Y_INIT
    for i in range(1, TS.size):
        t = (i - 1) * DT
        y = ys[i - 1]
        dw = dW(DT)
        del t
        # Sum up terms as in the Milstein method
        ys[i] = y + \
            Model.mu * y * DT + \
            Model.sigma * y * dw + \
            (Model.sigma**2 / 2) * y * (dw**2 - DT)

    return TS, ys


def plot_simulations(num_sims: int):
    """Plot several simulations in one image."""
    for _ in range(num_sims):
        plt.plot(*run_simulation())

    plt.xlabel("time (s)")
    plt.ylabel("y")
    plt.grid()
    plt.show()


if __name__ == "__main__":
    NUM_SIMS = 2
    plot_simulations(NUM_SIMS)
