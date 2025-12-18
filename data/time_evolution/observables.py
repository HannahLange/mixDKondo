import pyten as ptn
import os
import helper as H
import numpy as np

def energy(folder, state, lat):
    """
    Calculates the energy.
    """
    H.log("Calculate energy.")
    E = ptn.mp.expectation(state, lat.get("H_FH"))
    with open(folder+"E"+".csv", "w") as f:
        f.write(str(E))
    return E

def density(folder, length, state, lat):
    """
    Evaluates the local charge density.
    """
    H.log("Calculate density.")
    density = np.zeros((length))
    for x in range(length):
        n = ptn.mp.expectation(state, lat.get("n", x))
        density[x] = (np.real(n))
    np.save(folder+"density", density) 
    print("density:", density, density.sum())
    return density

def doublon_density(folder, length, state, lat):
    """
    Evaluates the local doublon density.
    """
    H.log("Calculate doublon density.")
    density = np.zeros((length))
    for x in range(length):
        n = ptn.mp.expectation(state, lat.get("nu", x)*lat.get("nd", x))
        density[x] = (np.real(n))
    np.save(folder+"doublon_density", density) 
    print("density:", density, density.sum())
    return density


def mag(folder, length, state, lat):
    """
    Evaluates the local magnetization
    """
    H.log("Calculate magnetization.")
    density = np.zeros((length))
    for x in range(length):
        n = ptn.mp.expectation(state, lat.get("sz", x))
        density[x] = (np.real(n))
    np.save(folder+"mag", density) 
    print("magnetization:", density, density.sum())
    return density

def spin_corrs(folder, length, state, lat):
    """
    Evaluates the spin-spin correlations
    """
    H.log("Calculate spin-spin correlations.")
    density = np.zeros((length, length))
    for x in range(length):
        for xprime in range(length):
            ss = lat.get("sz", x)*lat.get("sz", xprime)
            ss += 0.5*(lat.get("sp", x)*lat.get("sm", xprime)+lat.get("sm", x)*lat.get("sp", xprime))
            density[x,xprime] = (np.real(ptn.mp.expectation(state,ss)))
    np.save(folder+"spin_corrs", density) 
    print("spin corrs:", density)
    return density


