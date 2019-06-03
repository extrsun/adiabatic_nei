"""
Purpose: make a comparison plot: NEI v.s. CIE.

In this file, create the ionic fraction plot

Written by sw, Jun 02, 2019
"""

# Import used modules
import pyatomdb
import pickle, os
import numpy as np
import astropy.io.fits as pyfits
from astropy.io import ascii
import matplotlib as mpl
import matplotlib.pyplot as plt

#system parameters
rootpath = os.getcwd()+'/'

#Read the condition file
confile = rootpath+'adia.exp_phy.info'
conditions = ascii.read(confile)
ncondi = len(conditions)

Zlist = [1,2,6,7,8,10,12,14,16,18,20,26,28]
cie_spec = pickle.load(open(rootpath+'cie_case/tcspec_cie.pkl','rb'))
nei_spec = pickle.load(open(rootpath+'nei_case/tcspec_nei.pkl','rb'))
cie_spec_comp = pickle.load(open(rootpath+'cie_case/tcspec_comp_cie.pkl','rb'))
nei_spec_comp = pickle.load(open(rootpath+'nei_case/tcspec_comp_nei.pkl','rb'))

ebins = nei_spec['ebins']
nbins = len(cie_spec[0,:])
del_ebins = np.zeros(nbins,dtype=float)
del_ebins[0:(nbins-1)] = ebins[1:] - ebins[0:(nbins-1)]
del_ebins[nbins-1] = del_ebins[nbins-2]

nei_tspec = np.zeros([ncondi,nbins], dtype=float)
for Z in Zlist:
  nei_tspec += nei_spec[Z]

nei_tspec_comp = np.zeros([ncondi,nbins], dtype=float)
for Z in Zlist:
  nei_tspec_comp += nei_spec_comp[Z]

condi_index = [26, 100, 125, 200, 250, 283, 342]
for i in condi_index:
  plt.loglog(ebins, nei_tspec[i,:]/del_ebins*1e14, drawstyle='steps', \
    label='NEI_total', color='r')
  plt.loglog(ebins, cie_spec[i,:]/del_ebins*1e14, drawstyle='steps', \
    label='CIE_total', color='b')
  plt.loglog(ebins, nei_tspec_comp[i,:]/del_ebins*1e14, drawstyle='steps', \
    label='NEI_continuum', color='m')
  plt.loglog(ebins, cie_spec_comp[i,:]/del_ebins*1e14, drawstyle='steps', \
    label='CIE_continuum', color='c')
  plt.xlabel('Energy (keV)')
  plt.ylabel('10$^{-14}$ Cts s$^{-1}$ cm$^3$ keV$^{-1}$')
  plt.legend(loc=0)
  plt.xlim([0.1,2.0])
  plt.ylim([1e-3,1e2])
  radius=conditions[i]['R']/3.0856780e+18
  plt.savefig(rootpath+'figures/comp_spec/testmodel_comp_r%4.2fpc.png' % radius)
  plt.close()