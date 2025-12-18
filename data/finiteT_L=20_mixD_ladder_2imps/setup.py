import pyten as ptn
import os
import numpy as np
import helper as h
import observables as O
import argparse
import genHam as H
import time
from datetime import datetime
#import matplotlib.pyplot as plt


# Parse input. If an input parameter is not specified, the one from the config file is taken. 
parser = argparse.ArgumentParser()
parser.add_argument("-l",       "--length"          , type=int   , default = 10   , help="length of physical system (int, long side)")
parser.add_argument("-t1",       "--hopping1"         , type=float , default = 1.   , help="charge hopping in x-dir")
parser.add_argument("-J1",       "--J1"               , type=float , default = 0.5  , help="spin exchange along chain")
parser.add_argument("-t2",       "--hopping2"         , type=float , default = 0.   , help="charge hopping in x-dir")
parser.add_argument("-J2",       "--J2"               , type=float , default = 0.  , help="spin exchange along chain")
parser.add_argument("-Jperp",   "--Jperp"           , type=float , default = 0.5  , help="spin exchange with impurities")
parser.add_argument("-N",       "--N"               , type=int   , default = 5    , help="number of particles in the system")
parser.add_argument("-dmrg",    "--dmrg"            , type=int   , default = 1    , help="runs dmrg if ==1")
parser.add_argument("-imps",    "--imps"            , type=int, nargs='+', default = [3,7], help="list with locations of impurities")
parser.add_argument("-temp",    "--temp"            , type=int, default = 1, help="1: finite T, 0: only ground state search")

args = parser.parse_args()
length   =  args.length          # length of ladder
J1       =  args.J1              # spin exchange along chain
J2       =  args.J2              # spin exchange along chain
Jperp    =  args.Jperp           # spin exchange between impurity and chain
t1       =  args.hopping1        # charge hopping in x-dir
t2       =  args.hopping2        # charge hopping in x-dir
N        =  args.N               # number of particles in the system
imps     =  args.imps            # list with locations of impurities
run_dmrg     =  {0:False, 1:True}[args.dmrg]   # True: runs dmrg, False: loads converged mps if it exists


current_folder = os.getcwd().split("/")[-1]
folder = "t1="+str(t1)+"_t2="+str(t2)+"_J1="+str(J1)+"_J2="+str(J2)+"_Jperp="+str(Jperp)+"_N="+str(N)
ext, imps_str="_locs",""
for imp in imps:
    ext += "_"+str(imp)
    imps_str += " "+str(imp)

# --------- run dmrg ground state search --------------------------
if run_dmrg:
    print("python3 run_DMRG.py -l "+str(length)+" -t1 "+str(t1)+" -J1 "+str(J1)+" -t2 "+str(t2)+" -J2 "+str(J2)+" -Jperp "+str(Jperp)+" -N "+str(N)+" -imps"+imps_str+" -gs 1")
    os.system("python3 run_DMRG.py -l "+str(length)+" -t1 "+str(t1)+" -J1 "+str(J1)+" -t2 "+str(t2)+" -J2 "+str(J2)+" -Jperp "+str(Jperp)+" -N "+str(N)+" -imps"+imps_str+" -gs 1")

if args.temp and run_dmrg:
  if True:
    print("python3 run_DMRG.py -l "+str(length)+" -t1 "+str(t1)+" -J1 "+str(J1)+" -t2 "+str(t2)+" -J2 "+str(J2)+" -Jperp "+str(Jperp)+" -N "+str(N)+" -imps"+imps_str)
    os.system("python3 run_DMRG.py -l "+str(length)+" -t1 "+str(t1)+" -J1 "+str(J1)+" -t2 "+str(t2)+" -J2 "+str(J2)+" -Jperp "+str(Jperp)+" -N "+str(N)+" -imps"+imps_str)

  # ------ load ground state and lattice ----------
  state = ptn.mp.MPS("../converged_states/"+current_folder+"/Tinf_"+folder+ext+".mps")
  lat   = ptn.mp.Lattice("../converged_states/"+current_folder+"/mixD_"+folder+ext+".lat")

  print("----------Start time evolution---------")
  lat = H.gen_tJ(length, t1, t2, J1, J2, Jperp, imps, lat)

  def det_bond_dimension(state):
    bond_dimension = 0
    for i in range(state.size()):
        site_dimension = state[i].getReducedDims()[2]
        if site_dimension > bond_dimension:
            bond_dimension = site_dimension
    return bond_dimension

  maxt = 4
  dt = 0.05
  trunc = 1e-8
  N = int(maxt/dt)
  maxStates=3000

  energy, Ts, den, mag, ssimp, ssperp, szszimp, szszperp, ss, szsz = [],[],[],[],[],[],[],[],[],[]

  Psi = state.copy()
  t_conf       = ptn.krylov.Conf() # for krylov
  t_conf.dt    = -dt*1j
  t_conf.tbeg = 0j
  t_conf.tend = -dt*1j
  t_conf.errTolerance = trunc
  t_conf.threshold = 1e-8

  Hs = []
  Hs.append(lat.get("H"))
  tdvp = ptn.mp.krylov.Evolver_A(lat.get("H"), Psi, t_conf.copy())

  times = [[0.0, 0.0]]
  tstart = time.perf_counter()
  time_evolved = Psi.copy()
  setup_tdvp=True
  for i in range(0, N):
    m_i = det_bond_dimension(tdvp.get_ref_psi())
    if i<2 and i>0 or m_i<500:
        t_conf.tbeg = -i*dt*1j
        t_conf.tend = -(i+1)*dt*1j
        tdvp = ptn.mp.krylov.Evolver_A(lat.get("H"), time_evolved, t_conf.copy()) # for krylov
    elif setup_tdvp:
        t_conf       = ptn.tdvp.Conf()
        t_conf.dt    = -dt*1j
        t_conf.maxt  = -(maxt-i*dt)*1j
        t_conf.trunc = ptn.Truncation(trunc, maxStates=maxStates)
        t_conf.mode  = ptn.tdvp.Mode.TwoSite
        tdvp = ptn.mp.tdvp.PTDVP(time_evolved, Hs, t_conf.copy())
        setup_tdvp=False
    current_step = (i+1) * dt 
    tstep = time.perf_counter()

    print(datetime.now(), f"- Do sweep starting with bond dim m =", det_bond_dimension(tdvp.get_ref_psi()))
    try:
        tdvp.evolve_in_subspace() # for krylov
    except:
        tdvp.do_step()

    time_evolved = tdvp.get_ref_psi()
    time_evolved.normalise()
    times = np.append(times, [[current_step, time.perf_counter() - tstart ]], axis=0)

    E = ptn.mp.expectation(time_evolved, lat.get("H"))
    T = 1/(2*current_step)
    Ts.append(T)
    energy.append(E)
    #time_evolved.save("../converged_states/"+current_folder+"/gs_"+folder+ext+"T="+str(T)+".mps")    
    print(datetime.now(), f"- Step {current_step} finished, elapsed step time is {time.perf_counter() - tstep} secs, energy="+str(E)+" at T="+str(T))

    den.append(O.density(folder+"/", length, time_evolved, lat, ext, False))
    mag.append(O.mag(folder+"/", length+len(imps), time_evolved, lat, ext, False))
    ssimp.append(O.spin_corrs_imp(folder+"/", imps, time_evolved, lat, ext, False))
    ssperp.append(O.spin_corrs_perp(folder+"/", length,time_evolved, lat, ext, False))
    szszimp.append(O.sz_corrs_imp(folder+"/", imps, time_evolved, lat, ext, False))
    szszperp.append(O.sz_corrs_perp(folder+"/", length,time_evolved, lat, ext, False))
    szsz.append(O.sz_corrs_nn(folder+"/", length,time_evolved, lat, ext, False))
    ss.append(O.ss_corrs_nn(folder+"/", length,time_evolved, lat, ext, False))
    if i<10 or i%20==0:
      time_evolved.save("../converged_states/"+current_folder+"/gs_"+folder+ext+"T="+str(T)+".mps")
      np.save(folder+"/Es"+ext+".npy",energy)
      np.save(folder+"/Ts"+ext+".npy",Ts)
      np.save(folder+"/den"+ext+".npy",den)
      np.save(folder+"/mag"+ext+".npy",mag)
      np.save(folder+"/ssimp"+ext+".npy",ssimp)
      np.save(folder+"/ssperp"+ext+".npy",ssperp)
      np.save(folder+"/szszimp"+ext+".npy",szszimp)
      np.save(folder+"/szszperp"+ext+".npy",szszperp)
      np.save(folder+"/szsznn"+ext+".npy",szsz)
      np.save(folder+"/ssnn"+ext+".npy",ss)
  time_evolved.save("../converged_states/"+current_folder+"/gs_"+folder+ext+"T="+str(T)+".mps")
  np.save(folder+"/Es"+ext+".npy",energy)
  np.save(folder+"/Ts"+ext+".npy",Ts)
  np.save(folder+"/den"+ext+".npy",den)
  np.save(folder+"/mag"+ext+".npy",mag)
  np.save(folder+"/ssimp"+ext+".npy",ssimp)
  np.save(folder+"/ssperp"+ext+".npy",ssperp)
  np.save(folder+"/szszimp"+ext+".npy",szszimp)
  np.save(folder+"/szszperp"+ext+".npy",szszperp)
  np.save(folder+"/ssnn"+ext+".npy",ss)
else:
  files = os.listdir("../converged_states/"+current_folder)
  files = [file for file in files if "T=" in file and folder+ext in file]
  Ts = sorted([float(file.split("T=")[1].split(".mps")[0]) for file in files],reverse=True)
  Ts = [T for T in Ts if T>=1/(2*4)]
  print(Ts)
  ss, szsz = [],[]
  for T in Ts:
    print("T=",T)
    state = ptn.mp.MPS("../converged_states/"+current_folder+"/gs_"+folder+ext+"T="+str(T)+".mps")
    lat   = ptn.mp.Lattice("../converged_states/"+current_folder+"/mixD_"+folder+ext+".lat")
    ss.append(np.array([O.ss_corrs_nn(folder+"/", length,state, lat, ext, False),T]))
    szsz.append(np.array([O.sz_corrs_nn(folder+"/", length,state, lat, ext, False),T]))
  np.save(folder+"/ss_nn"+ext+".npy",np.array(ss))
  np.save(folder+"/szsz_nn"+ext+".npy",np.array(szsz))
