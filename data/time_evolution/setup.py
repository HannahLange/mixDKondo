import pyten as ptn
import os
import numpy as np
import helper as h
import observables as O
import argparse
import genHam as H
import time
from datetime import datetime
ptn.threading.setTensorNum(16)

# Parse input. If an input parameter is not specified, the one from the config file is taken. 
parser = argparse.ArgumentParser()
parser.add_argument("-l",       "--length"          , type=int   , default = 21   , help="length of physical system (int, long side)")
parser.add_argument("-t",       "--hopping"         , type=float , default = 1.   , help="charge hopping in x-dir")
parser.add_argument("-timp",    "--hopping_imp"     , type=float , default = 1.   , help="charge hopping to and from impurity site")
parser.add_argument("-U",       "--U"               , type=float , default = -2.5 , help="on-site repulsion U")
parser.add_argument("-Uimp",    "--Uimp"            , type=float , default = -7.5 , help="on-site repulsion U at impurity site")
parser.add_argument("-Delta",   "--Delta"           , type=float , default = 8.  , help="potential offset of impurity site")
parser.add_argument("-N",       "--N"               , type=int   , default = 21   , help="number of particles in the system")
parser.add_argument("-dmrg",    "--dmrg"            , type=int   , default = 1    , help="runs dmrg if ==1")

args = parser.parse_args()
length   =  args.length          # length of ladder
t        =  args.hopping         # hopping in x-dir
timp     =  args.hopping_imp     # hopping to and from impurity site
U        =  args.U               # on-site repulsive potential
Uimp     =  args.Uimp            # on-site repulsive potential at impurity site
Delta    =  args.Delta           # potential offset of impurity site
N        =  args.N               # number of particles in the system
run_dmrg     =  {0:False, 1:True}[args.dmrg]   # True: runs dmrg, False: loads converged mps if it exists


current_folder = os.getcwd().split("/")[-1]
folder = "U="+str(U)+"_Uimp="+str(Uimp)+"_Delta="+str(Delta)+"_t="+str(t)+"_timp="+str(timp)+"_N="+str(N)
folder0 = "init_U="+str(U)+"_t="+str(t)+"_N="+str(N)

# --------- run dmrg ground state search --------------------------
if run_dmrg:
    print("python3 run_DMRG.py -l "+str(length)+" -t "+str(t)+" -timp "+str(timp)+" -Delta "+str(Delta)+" -U "+str(U)+" -Uimp "+str(Uimp)+" -N "+str(N))
    os.system("python3 run_DMRG.py -l "+str(length)+" -t "+str(t)+" -timp "+str(timp)+" -Delta "+str(Delta)+" -U "+str(U)+" -Uimp "+str(Uimp)+" -N "+str(N))

# ------ load ground state and lattice ----------
state = ptn.mp.MPS("../converged_states/"+current_folder+"/gs_"+folder0+".mps")
lat   = ptn.mp.Lattice("../converged_states/"+current_folder+"/mixD_"+folder0+".lat")

if run_dmrg:
    # ------ calculate observables ------------------
    den = O.density(folder0+"/", length, state, lat)
    O.doublon_density(folder0+"/", length, state, lat)
    O.energy(folder0+"/", state, lat)
    O.mag(folder0+"/", length, state, lat)
    O.spin_corrs(folder0+"/", length, state, lat)



print("----------Start time evolution---------")

mid = 0
lat = H.gen_mixD_FH(length, t, timp, U, Uimp, Delta, False)
def det_bond_dimension(state):
    bond_dimension = 0
    for i in range(state.size()):
        site_dimension = state[i].getReducedDims()[2]
        if site_dimension > bond_dimension:
            bond_dimension = site_dimension
    return bond_dimension

maxt = 40
dt = 0.04
trunc=1e-8
N = int(maxt/dt)
maxStates=1500

A_jt_0 = np.zeros(shape=(3, N+1, length), dtype=np.complex128)
A_jt_1 = np.zeros(shape=(3, N+1, length), dtype=np.complex128)
A_jt_2 = np.zeros(shape=(3, N+1, length), dtype=np.complex128)
for lattice_site in range(length):
    A_jt_0[0, 0, lattice_site] = lattice_site
    A_jt_0[1, 0, lattice_site] = 0
    A_jt_0[2, 0, lattice_site] = ptn.mp.expectation(state, lat.get("sz", lattice_site), state) 
    A_jt_1[0, 0, lattice_site] = lattice_site
    A_jt_1[1, 0, lattice_site] = 0
    A_jt_1[2, 0, lattice_site] = ptn.mp.expectation(state, lat.get("n", lattice_site), state)

for lattice_site in range(length):
  ss =  lat.get("sz", mid)*lat.get("sz", lattice_site)
  ss += 0.5*(lat.get("sp", mid)*lat.get("sm", lattice_site)+lat.get("sm", mid)*lat.get("sp", lattice_site))
  A_jt_2[0, 0, lattice_site] = lattice_site
  A_jt_2[1, 0, lattice_site] = 0
  A_jt_2[2, 0, lattice_site] = ptn.mp.expectation(state, ss, state)

Psi = state.copy()
t_conf       = ptn.krylov.Conf() # for krylov
t_conf.dt    = dt
t_conf.tbeg = 0
t_conf.tend = dt
t_conf.errTolerance = trunc
t_conf.threshold = 1e-8

Hs = []
Hs.append(lat.get("H_FH_dyn"))
tdvp = ptn.mp.krylov.Evolver_A(lat.get("H_FH_dyn"), Psi, t_conf.copy())


times = [[0.0, 0.0]]
tstart = time.perf_counter()

time_evolved=state
for i in range(0, N):
    if False: #i>0:
        t_conf       = ptn.krylov.Conf() # for krylov
        t_conf.dt    = dt
        t_conf.errTolerance = trunc 
        t_conf.threshold = 1e-8
        t_conf.tbeg = i*dt
        t_conf.tend = (i+1)*dt
        tdvp = ptn.mp.krylov.Evolver_A(lat.get("H_FH_dyn"), time_evolved, t_conf) # for krylov
    if i==0:
        t_conf       = ptn.tdvp.Conf()
        t_conf.dt    = dt
        t_conf.maxt  = maxt
        t_conf.trunc = ptn.Truncation(trunc, maxStates=maxStates)
        t_conf.mode  = ptn.tdvp.Mode.GSE
        t_conf.gse_conf.adaptive = False
        #t_conf.mode  = ptn.tdvp.Mode.TwoSite
        tdvp = ptn.mp.tdvp.PTDVP(time_evolved, Hs, t_conf.copy())
    current_step = (i+1) * dt
    tstep = time.perf_counter()

    print(datetime.now(), f"- Do sweep starting with bond dim m =", det_bond_dimension(tdvp.get_ref_psi()))
    try:
        tdvp.evolve_in_subspace() # for krylov
    except:
        tdvp.do_step()
                
    time_evolved = tdvp.get_ref_psi()
    #time_evolved = tdvp.get_psi(False)
    times = np.append(times, [[current_step, time.perf_counter() - tstart ]], axis=0)

    for lattice_site in range(length):
        if lattice_site==mid: print(ptn.mp.expectation(time_evolved, lat.get("sz", lattice_site)), ptn.mp.expectation(time_evolved, lat.get("n", lattice_site)))
        A_jt_0[0, i+1, lattice_site] = lattice_site
        A_jt_0[1, i+1, lattice_site] = current_step
        A_jt_0[2, i+1, lattice_site] = ptn.mp.expectation(time_evolved, lat.get("sz", lattice_site), time_evolved)
        A_jt_1[0, i+1, lattice_site] = lattice_site
        A_jt_1[1, i+1, lattice_site] = current_step
        A_jt_1[2, i+1, lattice_site] = ptn.mp.expectation(time_evolved, lat.get("n", lattice_site), time_evolved)
    for lattice_site in range(length):
      ss =  lat.get("sz", mid)*lat.get("sz", lattice_site)
      ss += 0.5*(lat.get("sp", mid)*lat.get("sm", lattice_site)+lat.get("sm", mid)*lat.get("sp", lattice_site))
      print(ptn.mp.expectation(time_evolved, ss, time_evolved))
      A_jt_2[0, i+1, lattice_site] = lattice_site
      A_jt_2[1, i+1, lattice_site] = 0
      A_jt_2[2, i+1, lattice_site] = ptn.mp.expectation(time_evolved, ss, time_evolved)

    if current_step % 1==0:
      np.save(folder+"/Sz_evolution.npy", A_jt_0)
      np.save(folder+"/n_evolution.npy", A_jt_1)
      np.save(folder+"/Sigma_evolution.npy", A_jt_2)
    print(datetime.now(), f"- Step {current_step} finished, elapsed step time is {time.perf_counter() - tstep} secs")


np.save(folder+"/Sz_evolution.npy", A_jt_0)
np.save(folder+"/n_evolution.npy", A_jt_1)
np.save(folder+"/Sigma_evolution.npy", A_jt_2)
