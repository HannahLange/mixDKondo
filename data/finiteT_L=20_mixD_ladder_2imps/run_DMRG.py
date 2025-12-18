import pyten as ptn
import numpy as np
import os
import argparse

import helper as h
import genHam as H
import observables as O

# Parse input. If an input parameter is not specified, the one from the config file is taken. 
parser = argparse.ArgumentParser()
parser.add_argument("-l",       "--length"          , type=int   , default = 10   , help="length of physical system (int, long side)")
parser.add_argument("-t1",       "--hopping1"         , type=float , default = 1.   , help="charge hopping in x-dir")
parser.add_argument("-J1",       "--J1"               , type=float , default = 0.5  , help="spin exchange along chain")
parser.add_argument("-t2",       "--hopping2"         , type=float , default = 0.   , help="charge hopping in x-dir")
parser.add_argument("-J2",       "--J2"               , type=float , default = 0.  , help="spin exchange along chain")
parser.add_argument("-Jperp",   "--Jperp"           , type=float , default = 0.5  , help="spin exchange with impurities")
parser.add_argument("-N",       "--N"               , type=int   , default = 6    , help="number of particles in the system")
parser.add_argument("-imps",    "--imps"            , type=int, nargs='+', default = [3,7], help="list with locations of impurities")
parser.add_argument("-gs",    "--gs"            , type=int, default = 0, help="0: entangler, 1: ground state search")

args = parser.parse_args()
length   =  args.length          # length of ladder
J1       =  args.J1               # spin exchange along chain
J2       =  args.J2               # spin exchange along chain
Jperp    =  args.Jperp           # spin exchange between impurity and chain
t1       =  args.hopping1         # charge hopping in x-dir
t2       =  args.hopping2         # charge hopping in x-dir
N        =  args.N               # number of particles in the system
imps     =  args.imps            # list with locations of impurities
  
current_folder = os.getcwd().split("/")[-1]
folder = "t1="+str(t1)+"_t2="+str(t2)+"_J1="+str(J1)+"_J2="+str(J2)+"_Jperp="+str(Jperp)+"_N="+str(N)

ext = "_locs"
for imp in imps:
    ext += "_"+str(imp)
print("---------------" +folder)

#--------------------------DMRG ground state setup-----------------------------#
#-----------------------------------------------------------------#

# Generate lattice
lat = H.gen_entangler_tJ(length,t2,imps)
lat = H.gen_tJ(length, t1, t2, J1, J2, Jperp, imps, lat)

# Save lattice
for new_fol in ["../init/", "../converged_states/"]:
    if not os.path.exists(new_fol):
        os.mkdir(new_fol)
    if not os.path.exists(new_fol+current_folder):
        os.mkdir(new_fol+current_folder)
lat_save_str = "../init/"+current_folder+"/lat_"+folder+ext+".lat"
lat.save(lat_save_str)
# Generate initial state

h.log("Generate initial state: N = "+str(N))

if os.path.exists("../converged_states/"+current_folder+"/gs_"+folder+ext+".mps") and args.gs==1:
    h.log("Load existing state and run DMRG.")
    init = ptn.mp.MPS("../converged_states/"+current_folder+"/gs_"+folder+ext+".mps")
else:
    if args.gs==1:
      Sztot = 0.0
      print(str(float(N))+", 0.0, "+str(float(N))+", 0.0, "+str(Sztot))
      init = ptn.mp.generateCompleteState(lat, str(float(N))+", 0.0, "+str(float(N))+", 0.0, "+str(Sztot))
      for i,x in enumerate(imps):
        idx_p = h.xy_to_snake(x,1,4)
        idx_a = h.xy_to_snake(x,3,4)
        if i%2==0:
            init *= lat.get("chu", idx_p)
            init *= lat.get("chd", idx_a)
        else:
            init *= lat.get("chd", idx_p)
            init *= lat.get("chu", idx_a)
    else:
      print("Generate inital state for entangler")
      init = ptn.mp.generateCompleteState(lat, "0.0, 0.0, 0.0, 0.0, 0.0")
      dist = length//N
      for i in range(N):
        x = i*dist
        if i%2:
          x = length-x
        print(i,x)
        idx_p = h.xy_to_snake(x,0,4)
        idx_a = h.xy_to_snake(x,2,4)
        if i%2:
          init *= (lat.get("chu", idx_p)*lat.get("chd", idx_a))
        else:
          init *= (lat.get("chd", idx_p)*lat.get("chu", idx_a)) 
        init.normalise()
      for i,x in enumerate(imps):
        idx_p = h.xy_to_snake(x,1,4)
        idx_a = h.xy_to_snake(x,3,4)
        if i%2==0:
            init *= lat.get("chu", idx_p)
            init *= lat.get("chd", idx_a)
        else:
            init *= lat.get("chd", idx_p)
            init *= lat.get("chu", idx_a)
        init.normalise()
# Energy of initial stale
E_init = ptn.mp.expectation(init, lat.get("H_ent"))
h.log("Initial energy: " + str(E_init))
n = 0
for x in range(length):
    idx = h.xy_to_snake(x,1,4)
    n += (ptn.mp.expectation(init, lat.get("n", idx)))
    print("local density:", ptn.mp.expectation(init,lat.get("sz", idx)))
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
bond_dim     = [100,    256,    512,   1024]
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
if args.gs==0: pdmrg = ptn.mp.dmrg.PDMRG(init, [lat.get("H_ent")], dmrgconf)
else: pdmrg = ptn.mp.dmrg.PDMRG(init, [lat.get("H")], dmrgconf)

# run stages
E = []
for i in range(num_stages):
    mps_i = pdmrg.run()
    if args.gs==0: E.append(ptn.mp.expectation(mps_i, lat.get("H_ent")))
    else: E.append(ptn.mp.expectation(mps_i, lat.get("H")))
print("Final energy: ", ptn.mp.expectation(mps_i, lat.get("H_ent")), "(entangler), ",ptn.mp.expectation(mps_i, lat.get("H")), "(ground state)")
if args.gs==1: np.save(folder+"/E_gs"+ext+".npy", E[-1])
if not os.path.exists(folder):
    os.mkdir(folder)
with open(folder+"/info.txt", "a") as f:
    f.write("Delta E = "+str(np.real(E[-2]-E[-1]))+"\n")
    f.write("Sz_tot = "+str(sum([ptn.mp.expectation(mps_i, lat.get("sz",i)) for i in range(length)]))+"\n")
    f.write("SzSz-SxSx = "+str(sum([ptn.mp.expectation(mps_i, lat.get("sz",i)*lat.get("sz",i+1))-0.25*ptn.mp.expectation(mps_i, (lat.get("sp",i)+lat.get("sm",i))*(lat.get("sp",i+1)+lat.get("sm",i+1))) for i in range(length-1)]))+"\n")
# save state after last stage
if not os.path.exists("../converged_states/"+current_folder):
    os.mkdir("../converged_states/"+current_folder)
if args.gs==1:
  print("../converged_states/"+current_folder+"/gs_"+folder+ext+".mps")
  mps_i.save("../converged_states/"+current_folder+"/gs_"+folder+ext+".mps")
else:
  print("../converged_states/"+current_folder+"/Tinf_"+folder+ext+".mps")
  mps_i.save("../converged_states/"+current_folder+"/Tinf_"+folder+ext+".mps")
lat.save("../converged_states/"+current_folder+"/mixD_"+folder+ext+".lat")
for x in range(length):
    idx = h.xy_to_snake(x,1,4)
    n += (ptn.mp.expectation(mps_i, lat.get("n", idx)))
    print("local density:", ptn.mp.expectation(mps_i,lat.get("n", idx)))
print("density:", n)

if args.gs==1:
  ext += "_gs"
  den = O.density(folder+"/", length, mps_i, lat, ext)
  O.mag(folder+"/", length, mps_i, lat, ext)
  O.spin_corrs_imp(folder+"/", imps, mps_i, lat, ext)
  O.spin_corrs_perp(folder+"/", length, mps_i, lat, ext)
  O.sz_corrs_imp(folder+"/", imps, mps_i, lat, ext)
  O.sz_corrs_perp(folder+"/", length, mps_i, lat, ext)
else:
  O.check_infinite_T_state(folder+"/", length, mps_i, lat, ext, save=False)
