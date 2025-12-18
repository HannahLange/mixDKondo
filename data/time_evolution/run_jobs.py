import numpy as np
import os
import helper as h


run_dmrg = True
length   = 21        # length of system
t        = 1.0
U        = 0.
N        = 14+1

current_folder = os.getcwd()

dmrg = {True: 1, False: 0}[run_dmrg]
t = float(t)
U = float(U)
for Delta in [-4.]: #-20.0, 20.0, -12.0,-8.0,-4.0,4.0,8.0,12.0]:
  for Uimp in [8.0]: #-12.0,-8.0,-4.0,0.0, 4.0,8.0,12.0]:
    Delta = float(Delta)
    for timp in [1.0]:
      if True: #timp<np.abs(Delta+Uimp):
        J = (2*timp**2/(Uimp+Delta)+2*timp**2/(U-Delta))/t
        if True:
          timp = float(timp)
          folder = "U="+str(U)+"_Uimp="+str(Uimp)+"_Delta="+str(Delta)+"_t="+str(t)+"_timp="+str(timp)+"_N="+str(N)
          print(folder)
          print("J=",J)
          if not os.path.exists(folder):
              os.mkdir(folder) 
          file = "job.sh"
          with open(file, "r") as f:
            data = f.read()
          data = data.replace(".py", ".py -l "+str(length)+" -t "+str(t)+" -timp "+str(timp)+" -Delta "+str(Delta)+" -U "+str(U)+" -Uimp "+str(Uimp)+" -N "+str(N)+" -dmrg "+str(dmrg) )
          with open(folder+"/"+file, "w") as f:
            data = f.write(data)
          os.system("sbatch "+folder+"/"+file)
