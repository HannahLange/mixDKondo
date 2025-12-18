import pyten as ptn
import os
import helper as H
import numpy as np

def energy(folder, state, lat, ext):
    """
    Calculates the energy.
    """
    H.log("Calculate energy.")
    E = ptn.mp.expectation(state, lat.get("H"))
    with open(folder+"E"+ext+".csv", "w") as f:
        f.write(str(E))
    return E

def density(folder, length, state, lat, ext, save=True):
    """
    Evaluates the local charge density.
    """
    H.log("Calculate density.")
    density = np.zeros((2,length))

    for x in range(length):
        for y in range(2):
            idx = H.xy_to_snake(x,y,4)
            n = ptn.mp.expectation(state, lat.get("n", idx))
            density[y,x] = (np.real(n))
    if save: np.save(folder+"density"+ext, density) 
    print("density:", density, density.sum())
    return density


def mag(folder, length, state, lat,ext,save=True):
    """
    Evaluates the local magnetization
    """
    H.log("Calculate magnetization.")
    density = np.zeros((2,length))
    for x in range(length):
        for y in range(2):
            idx = H.xy_to_snake(x,y,4)
            n = ptn.mp.expectation(state, lat.get("sz", idx))
            density[y,x] = (np.real(n))
    if save: np.save(folder+"mag"+ext, density) 
    print("magnetization:", density, density.sum())
    return density

def spin_corrs_perp(folder, length, state, lat, ext, save=True):
    """
    Evaluates the spin-spin correlations
    """
    H.log("Calculate spin-spin correlations.")
    density = np.zeros((length))
    for x in range(length):
        idx, idxprime = H.xy_to_snake(x,0,4), H.xy_to_snake(x,1,4)
        ss = lat.get("sz", idx)*lat.get("sz", idxprime)
        ss += 0.5*(lat.get("sp", idx)*lat.get("sm", idxprime)+lat.get("sm", idx)*lat.get("sp", idxprime))
        density[x] = (np.real(ptn.mp.expectation(state,ss)))
    if save: np.save(folder+"spin_corrs_perp"+ext, density) 
    print("spin corrs perp:", density)
    return density

def sz_corrs_perp(folder, length, state, lat, ext, save=True):
    """
    Evaluates the spin-spin correlations
    """
    H.log("Calculate spin-spin correlations.")
    density = np.zeros((length))
    for x in range(length):
        idx, idxprime = H.xy_to_snake(x,0,4), H.xy_to_snake(x,1,4)
        ss = lat.get("sz", idx)*lat.get("sz", idxprime)
        density[x] = (np.real(ptn.mp.expectation(state,ss)))
    if save: np.save(folder+"sz_corrs_perp"+ext, density)
    print("sz corrs perp:", density)
    return density


def sz_corrs_nn(folder, length, state, lat, ext, save=True):
    """
    Evaluates the spin-spin correlations
    """
    H.log("Calculate spin-spin correlations.")
    density = np.zeros((length-1))
    for x in range(length-1):
        idx, idxprime = H.xy_to_snake(x,0,4), H.xy_to_snake(x+1,0,4)
        ss = lat.get("sz", idx)*lat.get("sz", idxprime)
        density[x] = (np.real(ptn.mp.expectation(state,ss)))
    if save: np.save(folder+"sz_corrs_nn"+ext, density)
    print("sz corrs NN:", density)
    return density

def ss_corrs_nn(folder, length, state, lat, ext, save=True):
    """
    Evaluates the spin-spin correlations
    """
    H.log("Calculate spin-spin correlations.")
    density = np.zeros((length-1))
    for x in range(length-1):
        idx, idxprime = H.xy_to_snake(x,0,4), H.xy_to_snake(x+1,0,4)
        ss = lat.get("sz", idx)*lat.get("sz", idxprime)
        ss += 0.5*(lat.get("sp", idx)*lat.get("sm", idxprime)+lat.get("sm", idx)*lat.get("sp", idxprime)) 
        density[x] = (np.real(ptn.mp.expectation(state,ss)))
    if save: np.save(folder+"ss_corrs_nn"+ext, density)
    print("ss corrs NN:", density)
    return density

def spin_corrs_imp(folder, locs, state, lat, ext, save=True):
    """
    Evaluates the spin-spin correlations
    """
    H.log("Calculate spin-spin correlations.")

    if len(locs)==2:
        x,xprime = H.xy_to_snake(locs[0],1,4), H.xy_to_snake(locs[1],1,4)
        ss = lat.get("sz", x)*lat.get("sz", xprime)
        ss += 0.5*(lat.get("sp", x)*lat.get("sm", xprime)+lat.get("sm", x)*lat.get("sp", xprime))
        density = (np.real(ptn.mp.expectation(state,ss)))
    if save: np.save(folder+"spin_corrs_imp"+ext, density)
    print("spin corrs imp:", density)
    return density

def sz_corrs_imp(folder, locs, state, lat, ext, save=True):
    """
    Evaluates the spin-spin correlations
    """
    H.log("Calculate spin-spin correlations.")

    if len(locs)==2:
        x,xprime = H.xy_to_snake(locs[0],1,4), H.xy_to_snake(locs[1],1,4)
        ss = lat.get("sz", x)*lat.get("sz", xprime)
        density = (np.real(ptn.mp.expectation(state,ss)))
    if save: np.save(folder+"sz_corrs_imp"+ext, density)
    print("sz corrs imp:", density)
    return density




def check_infinite_T_state(folder, length, state, lat, ext, save=True):
    """
    Evaluates the spin-spin correlations between physical sites and ancillas
    """
    H.log("Calculate spin-spin correlations between physical sites and ancillas.")
    density = np.zeros((2,length))
    for x in range(length):
        for y in range(2):
            idx, idxprime = H.xy_to_snake(x,y+2,4), H.xy_to_snake(x,y+2,4)
            ss = lat.get("sz", idx)*lat.get("sz", idxprime)
            ss += 0.5*(lat.get("sp", idx)*lat.get("sm", idxprime)+lat.get("sm", idx)*lat.get("sp", idxprime))
            density[y,x] = (np.real(ptn.mp.expectation(state,ss)))
    if save: np.save(folder+"spin_corrs_perp"+ext, density) 
    print("singlets?:", density)
    return density

