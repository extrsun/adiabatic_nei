"""
Purpose: plot some derived physical parameters under pre-defined
  grids: tau, accumulative tau, delta_r, etc.

Written by sw, Jun 19, 2019, initial release
"""

import pickle, os, math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# System parameters
rootpath=os.getcwd() + '/'
modname = ['300_log', '3000_log', '300_linear']
nmod = 3

# Read the physical grid
griddat = {}
for Z in range(0,nmod):
  gridread = pickle.load(open(rootpath+'adia_exp.phy_'+modname[Z]+'.pkl','rb'))
  griddat[Z] = gridread
# ngrid1 = len(griddat1)
# ngrid2 = len(griddat2)
# ngrid3 = len(griddat3)

# Deal with the dynamical timescale
arr_r = {}
for Z in range(0,nmod):
  arr_r[Z] = griddat[Z]['R']
del_r = {}
for Z in range(0,nmod):
  ngrid = len(arr_r[Z])
  delr_val = np.zeros(ngrid, dtype=float)
  delr_val[1:] = arr_r[Z][1:] - arr_r[Z][0:(ngrid-1)]
  del_r[Z] = delr_val
meanv_arr = {}
for Z in range(0,nmod):
  ngrid = len(arr_r[Z])
  meanv = np.zeros(ngrid, dtype=float)
  meanv[0]  = griddat[Z]['velo'][0]
  meanv[1:] = (griddat[Z]['velo'][1:] + griddat[Z]['velo'][0:(ngrid-1)])/2
  meanv_arr[Z] = meanv
meandens_arr = {}
for Z in range(0,nmod):
  ngrid = len(arr_r[Z])
  meandens = np.zeros(ngrid, dtype=float)
  meandens[0]  = griddat[Z]['dens'][0]
  meandens[1:] = (griddat[Z]['dens'][1:] + griddat[Z]['dens'][0:(ngrid-1)])/2
  meandens_arr[Z] = meandens

del_t = {}
for Z in range(0,nmod):
  del_t[Z] = del_r[Z]/meanv_arr[Z]
del_tau = {}
for Z in range(0,nmod):
  del_tau[Z] = del_t[Z] * meandens_arr[Z]
accum_tau = {}
for Z in range(0,nmod):
  accum_tau[Z] = np.cumsum(del_tau[Z])

for Z in range(0,nmod):
  arr_r[Z] /= 3.0856780e+18

mpl.style.use('default')
for Z in range(0,nmod):
  plt.plot(arr_r[Z][1:], del_r[Z][1:]/3.0856780e+18, color="C"+"%s"%(Z), \
    linewidth=1.5, label=modname[Z])
plt.xlabel('Radius (pc)')
plt.ylabel('$\Delta$r (pc)')
plt.ylim(4e-4,4e-1)
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=0)
plt.savefig(rootpath+'figures/complot_delta_r.eps')
plt.close()

for Z in range(0,nmod):
  plt.plot(arr_r[Z][1:], del_tau[Z][1:]/1e7, color="C"+"%s"%(Z), \
    linewidth=1.5, label=modname[Z])
plt.xlabel('Radius (pc)')
plt.ylabel('Travel Timescale (10$^{7}$ cm$^{-3}$ s)')
# plt.ylim(0.0,1.2)
plt.xscale('log')
# plt.yscale('log')
plt.legend(loc=0)
plt.savefig(rootpath+'figures/complot_tauval.eps')
plt.close()

for Z in range(0,nmod):
  plt.plot(arr_r[Z][1:], accum_tau[Z][1:]/1e8, color="C"+"%s"%(Z), \
    linewidth=1.5, label=modname[Z])
plt.xlabel('Radius (pc)')
plt.ylabel('Cumulative Travel Timescale (10$^{8}$ cm$^{-3}$ s)')
# plt.ylim(0.1,1e2)
plt.xscale('log')
plt.yscale('log')
plt.legend(loc=0)
plt.savefig(rootpath+'figures/complot_accum_tau.eps')
plt.close()

