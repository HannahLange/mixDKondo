import pyten as ptn
import numpy as np
import os
import argparse

import helper as h
import genHam as H

# Parse input. If an input parameter is not specified, the one from the config file is taken. 
parser = argparse.ArgumentParser()
parser.add_argument("-l",       "--length"          , type=int   , default = 21   , help="length of physical system (int, long side)")
parser.add_argument("-t",       "--hopping"         , type=float , default = 1.   , help="charge hopping in x-dir")
parser.add_argument("-timp",	"--hopping_imp"     , type=float , default = 1.   , help="charge hopping to and from impurity site")
parser.add_argument("-U",       "--U"               , type=float , default = -2.5 , help="on-site repulsion U")
parser.add_argument("-Uimp",	"--Uimp"            , type=float , default = -7.5 , help="on-site repulsion U at impurity site")
parser.add_argument("-Delta",   "--Delta"           , type=float , default = 0.   , help="potential offset of impurity site")
parser.add_argument("-N",       "--N"               , type=int   , default = 21   , help="number of particles in the system")

args = parser.parse_args()
length   =  args.length          # length of ladder
t        =  args.hopping         # hopping in x-dir
timp     =  args.hopping_imp     # hopping to and from impurity site
U        =  args.U               # on-siterepulsive potential
Uimp     =  args.Uimp            # on-siterepulsive potential at impurity
Delta    =  args.Delta           # potential offset of impurity site
N        =  args.N               # number of particles in the system

current_folder = os.getcwd().split("/")[-1]
folder = "U="+str(U)+"_Uimp="+str(Uimp)+"_Delta="+str(Delta)+"_t="+str(t)+"_timp="+str(timp)+"_N="+str(N)

print("---------------" +folder)


#--------------------------DMRG ground state setup-----------------------------#
#-----------------------------------------------------------------#

# Generate lattice
Delta0 = 0.0
folder0 = "init_U="+str(U)+"_t="+str(t)+"_N="+str(N)
lat = H.gen_mixD_FH(length, t, timp, U, Uimp, Delta0, True)
# Save lattice
for new_fol in ["../init/", "../converged_states/"]:
    if not os.path.exists(new_fol):
        os.mkdir(new_fol)
    if not os.path.exists(new_fol+current_folder):
        os.mkdir(new_fol+current_folder)
lat_save_str = "../init/"+current_folder+"/lat_"+folder+".lat"
lat.save(lat_save_str)
# Generate initial state

h.log("Generate initial state: N = "+str(N))
mid = 0
if False: #os.path.exists("../converged_states/"+current_folder+"/gs_"+folder0+".mps"):
    h.log("Load existing state and run DMRG.")
    init = ptn.mp.MPS("../converged_states/"+current_folder+"/gs_"+folder0+".mps")
else:
    Sztot = (N%2)/2
    print(N, Sztot)
    init = ptn.mp.generateSampledState(lat, str(0)+","+str(0))
    # create impurity
    print(mid)
    init *=  lat.get("chu", mid)
    # create other particles
    distance = (length-1)//N
    print(distance)
    for n in range(N-1):
        x = n+1
        if n%2==0: #
            init *=  lat.get("chu", x)
        else:
            init *=  lat.get("chd", x)

# Energy of initial stale
E_init = ptn.mp.expectation(init, lat.get("H_FH"))
h.log("Initial energy: " + str(E_init))
n = 0
for x in range(length):
    n += (ptn.mp.expectation(init, lat.get("n", x)))
    print("local density:", ptn.mp.expectation(init,lat.get("sz", x)))
print("density:", n)

#--------------------------Run DMRG-----------------------------#
#---------------------------------------------------------------#
t = []
import time

# DMRG config object
dmrgconf = ptn.dmrg.DMRGConfig()

# prefix to be used for log files
dmrgconf.prefix = '../data/out/dmrg_out'

# DMRG stages
bond_dim     = [100,    256,    512,   512+256]
max_sweeps   = [ 40,     40,     20,     10]
trunc_weight = [1e-14,  1e-14,  1e-14,  1e-14]
trunc_thres  = [1e-14,  1e-14,  1e-14,  1e-14]
min_E_diff   = [1e-12,  1e-12,  1e-12,  1e-12]
version      = ['DMRG3S','DMRG3S','DMRG3S','2DMRG']
num_stages   = len(bond_dim)

# add stages to DMRG config
for ii in range(num_stages):
    dmrgconf.stages += [ptn.dmrg.DMRGStage( "(m " + str(int(bond_dim[ii])) + " v " + str(version[ii]) + " x " + str(int(max_sweeps[ii])) + " t " + str(trunc_thres[ii]) + " w " + str(trunc_weight[ii]) + " d " + str(min_E_diff[ii]) + ")")]



# PDMRG management object. Initialised with our init state, the Hamiltonian and config object
pdmrg = ptn.mp.dmrg.PDMRG(init, [lat.get("H_FH")], dmrgconf)


# run stages
E = []
for i in range(num_stages):
    mps_i = pdmrg.run()
    E.append(ptn.mp.expectation(mps_i, lat.get("H_FH")))

if not os.path.exists(folder0):
    os.mkdir(folder0)
with open(folder0+"/info.txt", "a") as f:
    f.write("Delta E = "+str(np.real(E[-2]-E[-1]))+"\n")
    f.write("Sz_tot = "+str(sum([ptn.mp.expectation(mps_i, lat.get("sz",i)) for i in range(length)]))+"\n")
    f.write("SzSz-SxSx = "+str(sum([ptn.mp.expectation(mps_i, lat.get("sz",i)*lat.get("sz",i+1))-0.25*ptn.mp.expectation(mps_i, (lat.get("sp",i)+lat.get("sm",i))*(lat.get("sp",i+1)+lat.get("sm",i+1))) for i in range(length-1)]))+"\n")
# save state after last stage
if not os.path.exists("../converged_states/"+current_folder):
    os.mkdir("../converged_states/"+current_folder)
print("../converged_states/"+current_folder+"/gs_"+folder0+".mps")
mps_i.save("../converged_states/"+current_folder+"/gs_"+folder0+".mps")
lat.save("../converged_states/"+current_folder+"/mixD_"+folder0+".lat")

