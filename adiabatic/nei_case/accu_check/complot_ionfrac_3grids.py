"""
Purpose: make a comparison plot: NEI using left v.s. right
temperature values

In this file, create the ionic fraction plot

Written by sw, Jun 14, 2019
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
modname = ['300_log', '300_linear', '3000_log']
nmod = 3
Zlist = [6]

for M in range(0,nmod):
  confile = rootpath+'adia_exp.phy_'+modname[M]+'.pkl'
  conditions = pickle.load(open(confile, 'rb'))
  ncondi = len(conditions['kT'])

  #1.Compare the ionic fraction========================================
  nei1_ionfrac = pickle.load(open(rootpath+'tionfrac_C_'+modname[M]+'_l.pkl','rb'))
  nei2_ionfrac = pickle.load(open(rootpath+'tionfrac_C_'+modname[M]+'_r.pkl','rb'))
  nei3_ionfrac = pickle.load(open(rootpath+'tionfrac_C_'+modname[M]+'_m.pkl','rb'))

  # Only compare the ionic fraction of C, N, O, Ne, Mg, and Si because of
  # the color table limit
  # Zlist = [6,7,8,10,12,14,16,26]

  radius = conditions['R'][0:ncondi]/3.0856780e+18
  mpl.style.use('default')
  for Z in Zlist:
    elsymb = pyatomdb.atomic.Ztoelsymb(Z)
    begin_ion = Z-5
    ion_range = range(begin_ion,Z+1)
    for i in ion_range:
      plt.plot(radius, nei1_ionfrac[Z][i,0:ncondi], color='C'+"%s"%(i-begin_ion), \
        linewidth=2, label=elsymb+' '+pyatomdb.atomic.int2roman(i+1))
      plt.plot(radius, nei2_ionfrac[Z][i,0:ncondi], color='C'+"%s"%(i-begin_ion), \
        linewidth=2.5, linestyle='dashed')
    plt.plot(radius, np.zeros(ncondi), linewidth=1.5, label='Left')
    plt.plot(radius, np.zeros(ncondi), linewidth=2.0, \
      linestyle='dashed', label='Right')
    plt.xlabel('Radius (pc)')
    plt.ylabel('Ionic Fraction')
    plt.ylim(1e-7,1.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc=0)
    plt.savefig(rootpath+'figures/comp_ionfrac_'+modname[M]+'_%s.eps' %elsymb)
    plt.close()
  
    for i in ion_range:
      plt.plot(radius, np.abs(nei1_ionfrac[Z][i,0:ncondi] - \
        nei2_ionfrac[Z][i,0:ncondi]), color='C'+"%s"%(i-begin_ion), \
        linewidth=1.0, label=elsymb+' '+pyatomdb.atomic.int2roman(i+1), \
        linestyle='solid')
      plt.plot(radius, np.abs(nei2_ionfrac[Z][i,0:ncondi] - \
        nei3_ionfrac[Z][i,0:ncondi]), color='C'+"%s"%(i-begin_ion), \
        linewidth=1.5, linestyle='dashed')
    plt.plot(radius, np.zeros(ncondi), linewidth=1.5, label='Left-Right')
    plt.plot(radius, np.zeros(ncondi), linewidth=2.0, \
      linestyle='dashed', label='Right-Median')
    plt.xlabel('Radius (pc)')
    plt.ylabel('Ionic Fraction Difference')
    plt.ylim(1e-10,1e-3)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc=0)
    plt.savefig(rootpath+'figures/comp_iondiff_'+modname[M]+'_%s.eps' %elsymb)
    plt.close()
  
