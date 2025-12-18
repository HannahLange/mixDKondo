import numpy as np
import os
import helper as h


run_dmrg = True
length   = 20        # length of system
t1       = 1.0       # hopping in x-dir
J1       = 0.5
t2       = 0.0 
J2       = np.round(J1*(t2/t1)**2,4)

current_folder = os.getcwd()
temp = 1

dmrg = {True: 1, False: 0}[run_dmrg]
t1,t2 = float(t1), float(t2)
J1,J2 = float(J1), float(J2)

if run_dmrg:
  for Jperp in [3.0]:
    for N in [6,8,10,12,14,16,18,20]:
      for d in range(1,length-5):
        imps = " 3 "+str(3+d)
        folder = "t1="+str(t1)+"_t2="+str(t2)+"_J1="+str(J1)+"_J2="+str(J2)+"_Jperp="+str(Jperp)+"_N="+str(N)
        print(folder)
        if not os.path.exists(folder):
          os.mkdir(folder) 
        file = "job.sh"
        with open(file, "r") as f:
          data = f.read()
        data = data.replace(".py", ".py -l "+str(length)+" -t1 "+str(t1)+" -t2 "+str(t2)+" -J1 "+str(J1)+" -J2 "+str(J2)+" -Jperp "+str(Jperp)+" -N "+str(N)+" -dmrg "+str(dmrg)+" -temp "+str(temp) +" -imps"+imps)
        with open(folder+"/"+file, "w") as f:
          data = f.write(data)
        os.system("sbatch "+folder+"/"+file)
else:
  for Jperp in [0.1,0.5,1.5,3.0]:
    for N in [6,8,10,12,14,16,18,20]:
      for d in range(1,length-4):
        imps = " 2 "+str(2+d)
        folder = "t1="+str(t1)+"_t2="+str(t2)+"_J1="+str(J1)+"_J2="+str(J2)+"_Jperp="+str(Jperp)+"_N="+str(N)
        print(folder)
        if not os.path.exists(folder):
          os.mkdir(folder) 
        file = "job.sh"
        with open(file, "r") as f:
          data = f.read()
        if d==1: s = ".py -l "+str(length)+" -t1 "+str(t1)+" -t2 "+str(t2)+" -J1 "+str(J1)+" -J2 "+str(J2)+" -Jperp "+str(Jperp)+" -N "+str(N)+" -dmrg "+str(dmrg)+" -temp "+str(temp) +" -imps"+imps
        else: s+= "\npython3 setup.py -l "+str(length)+" -t1 "+str(t1)+" -t2 "+str(t2)+" -J1 "+str(J1)+" -J2 "+str(J2)+" -Jperp "+str(Jperp)+" -N "+str(N)+" -dmrg "+str(dmrg)+" -temp "+str(temp) +" -imps"+imps
      data = data.replace(".py", s)
      with open(folder+"/"+file, "w") as f:
        data = f.write(data)
      os.system("sbatch "+folder+"/"+file)
