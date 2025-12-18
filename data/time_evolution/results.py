#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 09:32:51 2022

@author: hannah.lange
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as m
from scipy.optimize import curve_fit
from matplotlib.colors import LogNorm, Normalize
from matplotlib.colorbar import ColorbarBase
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colors as mcolors
from scipy.optimize import curve_fit
from scipy.fft import fft,fftfreq


# Some parameters to make the plots look nice
params = {
    "text.usetex": True,
    "font.family": "serif",
    "legend.fontsize": 12,
    "figure.figsize": (19, 3.),
    "axes.labelsize": 14,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "lines.linewidth": 3,
    "lines.markeredgewidth": 1,
    "lines.markersize": 1,
    "lines.marker": "o",
    "patch.edgecolor": "black",
}
plt.rcParams.update(params)
plt.style.use("seaborn-deep")

print(os.environ["PATH"])
os.environ["PATH"] += os.pathsep + '/opt/local/bin'
print(os.getenv("PATH"))

cm = plt.get_cmap('tab10') 
cm2 = plt.get_cmap('coolwarm') 
cm3 = plt.get_cmap('Oranges') 





def load_data( Uimp, Delta, timp, N, obs="Sz"):
    folder = "U="+str(U)+"_Uimp="+str(Uimp)+"_Delta="+str(Delta)+"_t="+str(timp)+"_timp="+str(timp)+"_N="+str(N)+"/"+obs+"_evolution.npy"
    print(folder)
    data = np.load(folder)
    
    return data



length   = 3      
t        = 1.0 
U        = 8.0
N        = length



"""for timp in [1.]:
    fig, ax = plt.subplots(1,10)
    for Delta in [4.0]:
      for u,Uimp in enumerate([8.]):
          Delta = float(Delta)   
          if timp<np.abs(Delta+Uimp):
              J = 2*timp**2/(Uimp+Delta)+2*timp**2/(-Delta)
              if True:
                  data = load_data(Uimp,Delta,timp,length, "n").real
                  plt.plot(data[2,:,length//2])
                  data2 = load_data(Uimp,Delta,timp,length, "Sigma").real
                  plt.plot(data2[2,:,:].sum(axis=-1)-data2[2,:,length//2])
                  plt.show()"""

for timp in [1.,1.5, 2.0]:
    t=timp
    for Delta in [0.0, 4.0]:
      for u,Uimp in enumerate([8.]):
          Delta = float(Delta)   
          if timp<np.abs(Delta+Uimp):
              J = 2*timp**2/(Uimp+Delta)+2*timp**2/(U-Delta)
              if True:
                  data = load_data(Uimp,Delta,timp,length, "n").real
                  plt.plot(data[1,:,0],data[2,:,length//2], label="density")
                  data2 = load_data(Uimp,Delta,timp,length, "Sigma").real
                  plt.plot(data[1,:,0],data2[2,:,:].sum(axis=-1)-data2[2,:,length//2], label="Sigma")
                  plt.show()           