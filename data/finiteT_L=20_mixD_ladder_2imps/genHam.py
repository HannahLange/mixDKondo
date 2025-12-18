import pyten as ptn
import numpy as np
import os

import helper as hlp

#------------------------------------------------------
#                                                      |
#        Lattices for ground state calculations        |
#                                                      |
#------------------------------------------------------

def gen_tJ(length, t1, t2, J1,J2, Jperp, imp_locs, lat=None,**kwargs):

    """
    MixD t-J model
    Arguments:   length     --    length of the system
                 t          --    hopping along the x-direction
                 J          --    exchange along chain
                 Jperp      --    exchange with impurities
    Returns:     lat        --    U(1)^(width) x U(1) symmetry conserving lattice file
    """
    ########################Lattice###########################
    hlp.log("Initializing lattice")
    if lat==None: lat = ptn.mp.lat.u1u1.gentJLadderMixD(length, 4)

    hlp.log("Setting up t-J Hamiltonian with two impurities")
    N_lat = length
    h = []
    imp_locs = [l+i for i,l in enumerate(imp_locs)]
    ## Bonds along legs
    hlp.log("Working on hopping terms")
    for x in range(length-1):
        idx1, idx2 = hlp.xy_to_snake(x,0,4),hlp.xy_to_snake(x+1,0,4)
        h += [-t1 * (  lat.get("cd", idx1) * lat.get("chd", idx2)  + lat.get("cu", idx1) * lat.get("chu", idx2)  ) ]
        h += [-t1 * (  lat.get("cd", idx2) * lat.get("chd", idx1)  + lat.get("cu", idx2) * lat.get("chu", idx1)  ) ]
    for x in range(length-1):
        idx1, idx2 = hlp.xy_to_snake(x,1,4),hlp.xy_to_snake(x+1,1,4)        
        h += [-t2 * (  lat.get("cd", idx1) * lat.get("chd", idx2)  + lat.get("cu", idx1) * lat.get("chu", idx2)  ) ]
        h += [-t2 * (  lat.get("cd", idx2) * lat.get("chd", idx1)  + lat.get("cu", idx2) * lat.get("chu", idx1)  ) ]

    # spin-spin interactions
    hlp.log("Working on spin exchange terms")
    for x in range(length-1):
        idx1, idx2 = hlp.xy_to_snake(x,0,4),hlp.xy_to_snake(x+1,0,4) 
        h += [ J1 * (  lat.get("sz", idx1) * lat.get("sz", idx2) - 0.25* lat.get("n", idx1) * lat.get("n", idx2) ) ]
        h += [0.5 * J1 * (  lat.get("sp", idx1) * lat.get("sm", idx2) + lat.get("sm", idx1) * lat.get("sp", idx2) ) ]
    for x in range(length-1):
        idx1, idx2 = hlp.xy_to_snake(x,1,4),hlp.xy_to_snake(x+1,1,4) 
        h += [ J2 * (  lat.get("sz", idx1) * lat.get("sz", idx2) - 0.25* lat.get("n", idx1) * lat.get("n", idx2) ) ]
        h += [0.5 * J2 * (  lat.get("sp", idx1) * lat.get("sm", idx2) + lat.get("sm", idx1) * lat.get("sp", idx2) ) ]
    for x in range(length):
        idx1, idx2 = hlp.xy_to_snake(x,0,4),hlp.xy_to_snake(x,1,4)
        h += [ Jperp * (  lat.get("sz", idx1) * lat.get("sz", idx2) - 0.25* lat.get("n", idx1) * lat.get("n", idx2) ) ]
        h += [0.5 * Jperp * (  lat.get("sp", idx1) * lat.get("sm", idx2) + lat.get("sm", idx1) * lat.get("sp", idx2) ) ] 

    hlp.log("Adding and truncating")
    h = ptn.mp.addLog(h)
    lat.add("H", "t-J Hamiltonian with impurities", h)
    return lat


def gen_entangler_tJ(length, t2, imp_locs, **kwargs):
    """
    MixD t-J model
    Arguments:   length     --    length of the system
                 t2         --    hopping along the x-direction in second chain
    Returns:     lat        --    U(1)^(2) x U(1) symmetry conserving lattice file
    """

    hlp.log("Initializing lattice")
    lat = ptn.mp.lat.u1u1.gentJLadderMixD(length, 4)

    hlp.log("Setting up t-J entangler")

    hgs = []
    y1, y2 = 0,0
    for x1 in range(length):
        for x2 in range(length):
            idx1_p, idx2_p = hlp.xy_to_snake(x1,y1,4),hlp.xy_to_snake(x2,y2,4)
            idx1_a, idx2_a = hlp.xy_to_snake(x1,y1+2,4),hlp.xy_to_snake(x2,y2+2,4)
            print(x1,y1,idx1_p,idx1_a)
            if idx1_p!=idx2_p:
                Delta_dag = (lat.get("chd", idx1_a) * lat.get("chu", idx1_p) - lat.get("chu", idx1_a) * lat.get("chd", idx1_p)  )
                Delta = (lat.get("cu", idx2_p) * lat.get("cd", idx2_a) - lat.get("cd", idx2_p) * lat.get("cu", idx2_a)  )
                hgs += [-1/2*lat.get("cu", idx2_p) * lat.get("cd", idx2_a)*lat.get("chd", idx1_a) * lat.get("chu", idx1_p)]
                hgs += [1/2*( lat.get("cd", idx2_p) * lat.get("cu", idx2_a)*lat.get("chd", idx1_a) * lat.get("chu", idx1_p))]
                hgs += [1/2*(lat.get("cu", idx2_p) * lat.get("cd", idx2_a)*lat.get("chu", idx1_a) * lat.get("chd", idx1_p) )]
                hgs += [-1/2*lat.get("cd", idx2_p) * lat.get("cu", idx2_a) *lat.get("chu", idx1_a) * lat.get("chd", idx1_p)]
                
    y1, y2 = 1, 1
    for x1 in range(length):
        for x2 in range(length):
            idx1_p, idx2_p = hlp.xy_to_snake(x1,y1,4),hlp.xy_to_snake(x2,y2,4)
            idx1_a, idx2_a = hlp.xy_to_snake(x1,y1+2,4),hlp.xy_to_snake(x2,y2+2,4)
            if t2!=0:
                if idx1_p!=idx2_p:
                    hgs += [-1/2*lat.get("cu", idx2_p) * lat.get("cd", idx2_a)*lat.get("chd", idx1_a) * lat.get("chu", idx1_p)]
                    hgs += [1/2*( lat.get("cd", idx2_p) * lat.get("cu", idx2_a)*lat.get("chd", idx1_a) * lat.get("chu", idx1_p))]
                    hgs += [1/2*(lat.get("cu", idx2_p) * lat.get("cd", idx2_a)*lat.get("chu", idx1_a) * lat.get("chd", idx1_p) )]
                    hgs += [-1/2*lat.get("cd", idx2_p) * lat.get("cu", idx2_a) *lat.get("chu", idx1_a) * lat.get("chd", idx1_p)]
            else:
                if x2==0 and x1 in imp_locs:
                    print(x1,y1,idx1_p,idx1_a)
                    hgs += [-1*(lat.get("sm", idx1_p) * lat.get("sp", idx1_a)+lat.get("sm", idx1_a) * lat.get("sp", idx1_p))]
    hlp.log("Adding and truncating")
    hgs = ptn.mp.addLog(hgs)
    lat.add("H_ent", "t-J entangler for conserved spin and charge", hgs)
    return lat
