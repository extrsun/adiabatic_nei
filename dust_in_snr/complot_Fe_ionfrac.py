"""
Purpose: make a comparison plot: Zehao's code v.s. my approach

In this file, create the ionic fraction plot

Written by sw, Jun 20, 2019
"""

# Import used modules
import pyatomdb
import pickle, os
import numpy as np
import astropy.io.fits as pyfits
import matplotlib as mpl
import matplotlib.pyplot as plt

def complot_Fe_ionfrac_module(log=False):
  #system parameters
  rootpath = os.getcwd()+'/'
  kTlist  = np.array([1.,2.,4.,10.]) # in keV
  n_kT = len(kTlist)

  Z = 26
  elsymb = pyatomdb.atomic.Ztoelsymb(Z)

  if log:
    subdirec = 'log'
  else:
    subdirec = 'linear'
  modname = ['Zehao', 'Mine']
  ionfracfile_name = ['data_package/ionfrac_95_v3.0.9.pkl', \
    'my_approach/'+subdirec+'/dustsnr_ionfrac_Fe.pkl']
  nmod = 2

  ionfrac1 = pickle.load(open(rootpath + ionfracfile_name[0], 'rb'))
  ionfrac2 = pickle.load(open(rootpath + ionfracfile_name[1], 'rb'))

  begin_ion = 17
  ion_range = range(begin_ion+1,27)

  for M in range(0,n_kT):
    plt.plot(ionfrac1['taulist'], ionfrac1['ionfrac'][M,:,0], \
      linewidth=1.0, label='Neutral', color='C0')
    for iZ in ion_range:
      plt.plot(ionfrac1['taulist'], ionfrac1['ionfrac'][M,:,iZ+1], \
        color='C'+"%s"%(iZ-begin_ion), linewidth=1.0, \
        label=elsymb+' '+pyatomdb.atomic.int2roman(iZ+1))
    plt.plot(ionfrac2['taulist'], ionfrac2['ionfrac'][M,0,:], \
      linewidth=2.0, linestyle='dashed', color='C0')
    for iZ in ion_range:
      plt.plot(ionfrac2['taulist'], ionfrac2['ionfrac'][M,iZ+1,:], \
        color='C'+"%s"%(iZ-begin_ion), linewidth=2.0, linestyle='dashed')
    plt.plot(ionfrac1['taulist'], np.zeros(len(ionfrac1['taulist'])), \
      linewidth=1.0, label=modname[0])
    plt.plot(ionfrac2['taulist'], np.zeros(len(ionfrac2['taulist'])), \
      linewidth=2.0, linestyle='dashed',label=modname[1])
    plt.xlabel('Ionization Timescale (cm$^{-3}$ s)')
    plt.ylabel('Ionic Fraction')
    plt.xlim(1e10,1e13)
    plt.ylim(1e-6,1.0)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc=0)
    plt.savefig(rootpath+'figures/'+subdirec+'/comp_ionfrac_%04.1f.eps' % kTlist[M])
    plt.close()
   
if __name__ == '__main__':
  complot_Fe_ionfrac_module(log=True)
  complot_Fe_ionfrac_module(log=False)