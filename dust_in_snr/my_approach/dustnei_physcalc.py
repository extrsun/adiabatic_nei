"""
The physcalc module contains the routine that calculates the ionic
fraction and the NEI X-ray spectrum based on the physical conditions
in an ASCII file.

V0.1 - Wei Sun, May 22, 2019: initial release
V0.2 - Wei Sun, May 28, 2019: separate the CIE and NEI case
V0.3 - Wei Sun, Jun 01, 2019: standardize for github uploading
V0.4 - Wei Sun, Jun 03, 2019: add keyword "dolines" to calculate
       pure continuum spectrum
"""

try:
  import astropy.io.fits as pyfits
except ImportError:
  import pyfits

import pickle, os
import numpy as np
import pyatomdb
from astropy.io import ascii
import astropy

def calc_nei_ionfrac(Zlist, condifile=False, outfilename=False, \
  init_file=False, ktval=False, rootpath=False):

  """
  Calculate the ionic fraction based on the physical conditions in
  an ASCII file. The only input parameter is the index (array) of
  the elements.

  Parameters
  ----------
  Zlist: [int]
    list of element nuclear charges

  Keywords
  --------
  condifile: string or dictionary
    the ASCII file containing the physical condition array. can also
    pass in a structure read from the ASCII file;
  outfilename: str
    the name of output pickle file recording the ionic fraction.
    The name of the output file is adopted as following sequence:
      1. specified by <outfilename>;
      2. adopted from <init_file>, if applicable;
      3. "tionfrac_Element.List.pkl".
  init_file: str or dictionary of ionic fraction
    the pickle file containing the ionic fraction at prior condition
    position. Could be the dictionary of loaded from the pickle file;

  Returns
  -------
  No return, but the pickle file is created/updated with derived ionic
  fraction at condition positions.

  """

  # System parameters
  atomdbpath = os.environ['ATOMDB']
  ionbalfile = atomdbpath+'APED/ionbal/v3.0.7_ionbal.fits'
  if not pyatomdb.util.keyword_check(rootpath):
    rootpath = os.getcwd()+'/'

  # Parameters related to the element list
  NZlist = len(Zlist)
  Zmax = np.max(Zlist)

  # Check the setting of the condition array
  if pyatomdb.util.keyword_check(condifile):
    # If it is a string, look for the file name and read it if exists
    if isinstance(condifile, str):
      confile = os.path.expandvars(rootpath+condifile)
      if not os.path.isfile(confile):
        print("*** ERROR: no such condition file %s. Exiting ***" \
          %(confile))
        return -1
      if confile.split('.')[-1] == 'pkl':
        conditions = pickle.load(open(confile,'rb'))
      else:
        conditions = ascii.read(confile)
    elif isinstance(condifile, astropy.table.table.Table):
      conditions = condifile
    elif isinstance(condifile, dict):
      conditions = condifile
    else:
      print("Unknown data type for condition file. Please pass a " \
        "string or an ASCIIList")
      return -1
  ncondi = len(conditions['kT'])
  if pyatomdb.util.keyword_check(ktval):
    conditions['kT'] = np.zeros(ncondi, dtype=float) + ktval

  # The final result - ionfrac
  ionfrac = {}
  for Z in Zlist:
    ionfrac[Z] = np.zeros([Z+1,ncondi],dtype=float)
  # Alternative way:
  #   ionfrac = [np.zeros([Z+1,ncondi],dtype=float) for Z in Zlist]
  # which does not have the "Z" index.

  # settings of the initial ionic fraction file
  if pyatomdb.util.keyword_check(init_file):
    # If it is a string, look for the file name and read it if exists
    if isinstance(init_file, str):
      initfile = os.path.expandvars(rootpath+init_file)
      if not os.path.isfile(initfile):
        print("*** ERROR: no such initial ionic fraction file %s. " \
          "Exiting ***" %(initfile))
        return -1
      add_ionfrac = pickle.load(open(init_file,'rb'))
    elif isinstance(init_file, dict):
      add_ionfrac = init_file
    else:
      print("Unknown data type for condition file. Please pass a " \
        "string or an DICTList")
      return -1

  # settings of the name of the output pickle file
  if not pyatomdb.util.keyword_check(outfilename):
    if pyatomdb.util.keyword_check(init_file) and \
         isinstance(init_file, str):
      outfilename = init_file
    else:
      outfilename = 'tionfrac_'
      for Z in Zlist:
        outfilename += pyatomdb.atomic.Ztoelsymb(Z)
      outfilename += '.pkl'

  # Initial ionic fraction: specified by init_file and begin_index,
  # or at r=0
  ion_init = {}
  for Z in Zlist:
    ion_init[Z] = add_ionfrac[Z][:,0] / add_ionfrac[Z][0,0]
    ionfrac[Z][:,0] = ion_init[Z]

  # Some tabled physical parameters
  temp_arr = conditions['kT']
  dens_arr = conditions['dens']
  radi_arr = conditions['R']
  velo_arr = conditions['velo']
  # Some derived physical parameters
  delr_arr     = np.zeros(ncondi, dtype=float)
  delr_arr[1:] = radi_arr[1:]-radi_arr[0:(ncondi-1)]
  meankt_arr     = np.zeros(ncondi, dtype=float)
  meankt_arr[0]  = temp_arr[0]
  meankt_arr[1:] = (temp_arr[1:]+temp_arr[0:(ncondi-1)])/2
  meanv_arr     = np.zeros(ncondi, dtype=float)
  meanv_arr[0]  = velo_arr[0]
  meanv_arr[1:] = (velo_arr[1:]+velo_arr[0:(ncondi-1)])/2
  time_arr = delr_arr / meanv_arr
  meandens     = np.zeros(ncondi, dtype=float)
  meandens[0]  = dens_arr[0]
  meandens[1:] = (dens_arr[1:]+dens_arr[0:(ncondi-1)])/2
  tau_arr  = time_arr * meandens
  
  # The radius/zone cycle
  condi_index = range(1,ncondi)
  for l in condi_index:
    print('For Zone-%03d: R=%10.3e:...' % (l, radi_arr[l]), \
      end='', flush=True)    
    # Calculate the ionic fraction (MOST IMPORTANT!)
    ionbal = {}
    for Z in Zlist:
      ionbal[Z] = pyatomdb.apec.solve_ionbal_eigen(Z, meankt_arr[l], \
                      init_pop=ion_init[Z], tau=tau_arr[l], teunit='keV')
      if pyatomdb.util.keyword_check(init_file):
        ionbal[Z] = add_ionfrac[Z][0,l-1] * ionbal[Z] + \
                    add_ionfrac[Z][:,l] - add_ionfrac[Z][:,l-1]
        ionbal[Z] /= np.sum(ionbal[Z])
        # print(add_ionfrac[Z][0,l],meankt_arr[l],tau_arr[l],ionbal[Z])
      ionfrac[Z][:,l] = ionbal[Z]
    ion_init = ionbal
    print('  finished.')

  # Save calculated ionic fraction as pickle file
  tmp = open(rootpath + outfilename, 'wb')
  pickle.dump(ionfrac,tmp)
  tmp.close()
  return 0