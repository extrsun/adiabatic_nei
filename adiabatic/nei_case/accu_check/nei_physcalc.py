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
  init_file=False, begin_index=False, end_index=False, rootpath=False):

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
  init_file: str or dictionary of ionic fraction
    the pickle file containing the ionic fraction at prior condition
    position. Could be the dictionary of loaded from the pickle file;
  begin_index: int
    the beginning index of the condition position, where the ionic
    fraction will be calculated outwards till <end_index>, def: 0;
  end_index: int
    the ending index of the condition position, where the calculation
    of the ionic fraction will be stopped, def: len(conditions);
  outfilename: str
    the name of output pickle file recording the ionic fraction.
    The name of the output file is adopted as following sequence:
      1. specified by <outfilename>;
      2. adopted from <init_file>, if applicable;
      3. "tionfrac_Element.List.pkl".

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
        conditions = pickle.load(open(confile, 'rb'))
      else:
        conditions = ascii.read(confile)
    elif isinstance(condifile, astropy.table.table.Table):
      conditions = condifile
    else:
      print("Unknown data type for condition file. Please pass a " \
        "string or an ASCIIList")
      return -1
  ncondi = len(conditions['kT'])

  # The final result - ionfrac
  ionfrac_l = {}
  ionfrac_r = {}
  ionfrac_m = {}
  for Z in Zlist:
    ionfrac_l[Z] = np.zeros([Z+1,ncondi],dtype=float)
    ionfrac_r[Z] = np.zeros([Z+1,ncondi],dtype=float)
    ionfrac_m[Z] = np.zeros([Z+1,ncondi],dtype=float)
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
      old_ionfrac = pickle.load(open(init_file,'rb'))
    elif isinstance(init_file, dict):
      old_ionfrac = init_file
    else:
      print("Unknown data type for condition file. Please pass a " \
        "string or an DICTList")
      return -1
    # Deal with the length of read ionic fraction
    for Z in Zlist:
      ncondi_old_ionfrac = len(old_ionfrac[Z][0,:])
      if ncondi_old_ionfrac < ncondi:
        ionfrac_l[Z][:,0:ncondi_old_ionfrac] = old_ionfrac[Z]
      else:
        ionfrac_l[Z] = old_ionfrac[Z][:,0:ncondi]

  # settings of the name of the output pickle file
  if not pyatomdb.util.keyword_check(outfilename):
    if pyatomdb.util.keyword_check(init_file) and \
         isinstance(init_file, str):
      outfilename = init_file
    else:
      outfilename = 'tionfrac_'
      for Z in Zlist:
        outfilename += pyatomdb.atomic.Ztoelsymb(Z)
      # outfilename += '.pkl'

  # Initial ionic fraction: specified by init_file and begin_index,
  # or at r=0
  ion_init_l = {}
  if pyatomdb.util.keyword_check(begin_index) and begin_index > 0 \
     and begin_index < ncondi:
    for Z in Zlist:
      ion_init_l[Z] = ionfrac[Z][:,begin_index]
  else:
    begin_index = 0
    te_init = conditions['kT'][begin_index]/pyatomdb.const.KBOLTZ #in K
    for Z in Zlist:
      ion_init_l[Z] = pyatomdb.atomdb.get_ionfrac(ionbalfile,
                      Z, te_init)
      ionfrac_l[Z][:,0] = ion_init_l[Z]
      ionfrac_r[Z][:,0] = ion_init_l[Z]
      ionfrac_m[Z][:,0] = ion_init_l[Z]
  ion_init_r = ion_init_l
  ion_init_m = ion_init_l

  # Deal with the ending index
  if (not pyatomdb.util.keyword_check(end_index)) or end_index < 0 \
     or end_index >= ncondi:
    end_index = ncondi-1

  condi_index = range(begin_index+1,end_index+1)

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
  for l in condi_index:
    print('For Zone-%03d: R=%10.3e:...' % (l, radi_arr[l]), \
      end='', flush=True)
    # if l < 1000:
    #   trans_kt = temp_arr[l-1]
    # else:
    #   trans_kt = temp_arr[l]

    # Calculate the ionic fraction (MOST IMPORTANT!)
    # ionbal = pyatomdb.apec.calc_full_ionbal(temp_arr[l], tau=tau_arr[l],
    #            init_pop=ion_init, Zlist=Zlist, teunit='keV', cie=False)

    # A test of ionic fraction difference using the electron temperature
    # at the beginning, middle, and ending points
    ionbal_m = {}
    for Z in Zlist:
      ionbal_m[Z] = pyatomdb.apec.solve_ionbal_eigen(Z, meankt_arr[l], \
                      init_pop=ion_init_m[Z], tau=tau_arr[l], teunit='keV')
      ionfrac_m[Z][:,l] = ionbal_m[Z]
    ion_init_m = ionbal_m
    ionbal_l = {}
    for Z in Zlist:
      ionbal_l[Z] = pyatomdb.apec.solve_ionbal_eigen(Z, temp_arr[l-1], \
                    init_pop=ion_init_l[Z], tau=tau_arr[l], teunit='keV')
      ionfrac_l[Z][:,l] = ionbal_l[Z]
    ion_init_l = ionbal_l
    ionbal_r = {}
    for Z in Zlist:
      ionbal_r[Z] = pyatomdb.apec.solve_ionbal_eigen(Z, temp_arr[l], \
                      init_pop=ion_init_r[Z], tau=tau_arr[l], teunit='keV')
      ionfrac_r[Z][:,l] = ionbal_r[Z]
    ion_init_r = ionbal_r
    # maxdiff_abs = 0.0
    # maxdiff_rel = 0.0
    # for Z in Zlist:
    #   maxdiff_abs = max([max(np.abs(ionbal_r[Z]-ionbal_l[Z])), maxdiff_abs])
    #   max_index = np.argmax(np.abs(ionbal_r[Z]-ionbal_l[Z]))
    #   maxrelative = np.abs(ionbal_r[Z][max_index]-ionbal_l[Z][max_index]) / \
    #                   max([ionbal_r[Z][max_index],ionbal_l[Z][max_index]])
    #   for iZ in range(0,Z+1):
    #     maxdiff_rel = max([np.abs(ionbal_r[Z][iZ]-ionbal_l[Z][iZ]) / \
    #                          max([ionbal_r[Z][iZ],ionbal_l[Z][iZ],1e-7]), \
    #                        maxdiff_rel])
    # print("%10.4e,%10.4e,%10.4e" % (maxdiff_abs,ionbal_l[Z][max_index],maxdiff_rel), \
    #       end='', flush=True)
    print('  finished.')

  # Save calculated ionic fraction as pickle file
  tmp = open(outfilename+'_l.pkl','wb')
  pickle.dump(ionfrac_l,tmp)
  tmp.close()
  tmp = open(outfilename+'_r.pkl','wb')
  pickle.dump(ionfrac_r,tmp)
  tmp.close()
  tmp = open(outfilename+'_m.pkl','wb')
  pickle.dump(ionfrac_m,tmp)
  tmp.close()
