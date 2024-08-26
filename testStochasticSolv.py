import numpy as np
import matplotlib.pyplot as plt
import brainpy as bp
sigma = 10; beta = 8/3; 
rho = 28;   p = 0.1

def lorenz_g(x, y, z, t):
    return p * x, p * y, p * z

def lorenz_f(x, y, z, t):
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return dx, dy, dz

lorenz = bp.sdeint(f=lorenz_f, 
                   g=lorenz_g, 
                   intg_type=bp.integrators.ITO_SDE,
                   wiener_type=bp.integrators.SCALAR_WIENER)
