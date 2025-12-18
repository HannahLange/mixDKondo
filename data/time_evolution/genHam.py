import pyten as ptn
import numpy as np
import os

import helper as hlp

#------------------------------------------------------
#                                                      |
#        Lattices for ground state calculations        |
#                                                      |
#------------------------------------------------------

def gen_mixD_FH(length, t, timp, U, Uimp, Delta, isolate, **kwargs):

    """
    MixD t-J model
    Arguments:   length     --    length of the system
                 t          --    hopping along the x-direction
                 U          --    on-site repulsion
                 Delta      --    potential offset of impurity
    Returns:     lat        --    U(1)^(width) x U(1) symmetry conserving lattice file
    """
    #assert length%2==1
    ########################Lattice###########################
    hlp.log("Initializing lattice")
    lat = ptn.mp.lat.u1u1.genFermiHubbard(length)

    hlp.log("Setting up mixD t-J Hamiltonian")
    N_lat = length
    h = []

    ## Bonds along legs
    mid = 0
    hlp.log("Working on hopping terms")
    for x in range(length-1):
        idx1, idx2 = x, x+1
        if (mid not in [idx1,idx2]):
            print(idx1, idx2)
            h += [-t * (  lat.get("cd", idx1) * lat.get("chd", idx2)  + lat.get("cu", idx1) * lat.get("chu", idx2)  ) ]
            h += [-t * (  lat.get("cd", idx2) * lat.get("chd", idx1)  + lat.get("cu", idx2) * lat.get("chu", idx1)  ) ]
        elif (not isolate):
            print("imp hopping:",idx1, idx2)
            h += [-timp * (  lat.get("cd", idx1) * lat.get("chd", idx2)  + lat.get("cu", idx1) * lat.get("chu", idx2)  ) ]
            h += [-timp * (  lat.get("cd", idx2) * lat.get("chd", idx1)  + lat.get("cu", idx2) * lat.get("chu", idx1)  ) ]
    # on-site repulsion
    hlp.log("Working on on-site terms")
    for x in range(length):
      if (mid!=x):
        h += [ U * (lat.get("nu", x) * lat.get("nd", x))]
      else:
        print("impurity site:", x)
        h += [ Uimp * (lat.get("nu", x) * lat.get("nd", x))]

    hlp.log("Working on impurity site")
    h += [Delta*(lat.get("nu", mid) + lat.get("nd", mid))]

    hlp.log("Adding and truncating")
    h = ptn.mp.addLog(h)
    if isolate:
        lat.add("H_FH", "MixD t-J Hamiltonian", h)
    else:
        lat.add("H_FH_dyn", "MixD t-J Hamiltonian for dynamics", h)
    return lat


def gen_mixD_tJ(length, t, timp, J, Jimp, Delta, isolate, **kwargs):

    """
    MixD t-J model
    Arguments:   length     --    length of the system
                 t          --    hopping along the x-direction
                 U          --    on-site repulsion
                 Delta      --    potential offset of impurity
    Returns:     lat        --    U(1)^(width) x U(1) symmetry conserving lattice file
    """
    #assert length%2==1
    ########################Lattice###########################
    hlp.log("Initializing lattice")
    lat = ptn.mp.lat.u1u1.gentJ(length)

    hlp.log("Setting up mixD t-J Hamiltonian")
    N_lat = length
    h = []

    ## Bonds along legs
    mid = length//2
    hlp.log("Working on hopping terms")
    for x in range(length-1):
        idx1, idx2 = x, x+1
        if (mid not in [idx1,idx2]):
            print(idx1, idx2)
            h += [-t * (  lat.get("cd", idx1) * lat.get("chd", idx2)  + lat.get("cu", idx1) * lat.get("chu", idx2)  ) ]
            h += [-t * (  lat.get("cd", idx2) * lat.get("chd", idx1)  + lat.get("cu", idx2) * lat.get("chu", idx1)  ) ]
            h += [J  * (  lat.get("sz", idx1) * lat.get("sz", idx2) ) ]
            h += [0.5*J * (  lat.get("sp", idx1) * lat.get("sm", idx2) + lat.get("sm", idx1) * lat.get("sp", idx2) ) ]
        elif (not isolate):
            print("imp hopping:",idx1, idx2)
            h += [-timp * (  lat.get("cd", idx1) * lat.get("chd", idx2)  + lat.get("cu", idx1) * lat.get("chu", idx2)  ) ]
            h += [-timp * (  lat.get("cd", idx2) * lat.get("chd", idx1)  + lat.get("cu", idx2) * lat.get("chu", idx1)  ) ]
            h += [Jimp  * (  lat.get("sz", idx1) * lat.get("sz", idx2) ) ]
            h += [0.5*Jimp * (  lat.get("sp", idx1) * lat.get("sm", idx2) + lat.get("sm", idx1) * lat.get("sp", idx2) ) ]
    # on-site repulsion

    hlp.log("Adding and truncating")
    h = ptn.mp.addLog(h)
    if isolate:
        lat.add("H_FH", "MixD t-J Hamiltonian", h)
    else:
        lat.add("H_FH_dyn", "MixD t-J Hamiltonian for dynamics", h)
    return lat
