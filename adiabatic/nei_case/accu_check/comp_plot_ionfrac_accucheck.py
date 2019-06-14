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
# Zlist = [1,2,6,7,8,10,12,14,16,18,20,26,28]

#Read the condition file
confile = rootpath+'../adia.exp_phy.info'
conditions1 = ascii.read(confile)
ncondi1 = len(conditions1)
conditions2 = ascii.read(rootpath+'../adia.exp_phy.info')
ncondi2 = len(conditions2)

#1.Compare the ionic fraction========================================
nei1_ionfrac = pickle.load(open(rootpath+'../tionfrac_nei.pkl','rb'))
nei2_ionfrac = pickle.load(open(rootpath+'../backup/tionfrac_nei.pkl','rb'))

# Only compare the ionic fraction of C, N, O, Ne, Mg, and Si because of
# the color table limit
# Zlist = [6,7,8,10,12,14,16,26]
Zlist = [6]

radius1 = conditions1[0:ncondi1]['R']/3.0856780e+18
radius2 = conditions2[0:ncondi2]['R']/3.0856780e+18
mpl.style.use('default')
for Z in Zlist:
  elsymb = pyatomdb.atomic.Ztoelsymb(Z)
  if Z == 26:
    begin_ion = 16
    ion_range = [16,17,18,19,20,21,22]
  else:
    begin_ion = Z-5
    ion_range = range(begin_ion,Z+1)
  for i in ion_range:
    plt.plot(radius1, nei1_ionfrac[Z][i,0:ncondi1], color='C'+"%s"%(i-begin_ion), \
      linewidth=2, label=elsymb+' '+pyatomdb.atomic.int2roman(i+1))
    plt.plot(radius2, nei2_ionfrac[Z][i,0:ncondi2], color='C'+"%s"%(i-begin_ion), \
      linewidth=2.5, linestyle='dashed')
  plt.plot(radius1, np.zeros(ncondi1), linewidth=2.0, label='updated')
  plt.plot(radius2, np.zeros(ncondi2), linewidth=2.5, \
    linestyle='dashed', label='original')
  plt.xlabel('Radius (pc)')
  plt.ylabel('Ionic Fraction')
  plt.ylim(1e-7,1.5)
  plt.xscale('log')
  plt.yscale('log')
  plt.legend(loc=0)
  plt.savefig(rootpath+'figures/comp_ionfrac_comp_%s.png' %elsymb)
  plt.close()
