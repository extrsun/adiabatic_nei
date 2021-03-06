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
import astropy
from astropy.io import ascii
from scipy.signal import find_peaks

def calc_nei_ionfrac(Zlist, condifile=False, diff_thres=False, \
  init_file=False, begin_index=False, end_index=False, outfilename=False, \
  rootpath=False):

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
  
  if not pyatomdb.util.keyword_check(diff_thres):
    diff_thres = 1e-7

  # Check the setting of the condition array
  if pyatomdb.util.keyword_check(condifile):
    # If it is a string, look for the file name and read it if exists
    if isinstance(condifile, str):
      confile = os.path.expandvars(rootpath+condifile)
      if not os.path.isfile(confile):
        print("*** ERROR: no such condition file %s. Exiting ***" \
          %(confile))
        return -1
      conditions = ascii.read(confile)
    elif isinstance(condifile, astropy.table.table.Table):
      conditions = condifile
    else:
      print("Unknown data type for condition file. Please pass a " \
        "string or an ASCIIList")
      return -1
  ncondi = len(conditions)

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
        ionfrac[Z][:,0:ncondi_old_ionfrac] = old_ionfrac[Z]
      else:
        ionfrac[Z] = old_ionfrac[Z][:,0:ncondi]

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
  if pyatomdb.util.keyword_check(begin_index) and begin_index >0 \
     and begin_index < ncondi:
    for Z in Zlist:
      ion_init[Z] = ionfrac[Z][:,begin_index]
  else:
    begin_index = 0
    te_init = conditions[begin_index]['kT']/pyatomdb.const.KBOLTZ #in K
    for Z in Zlist:
      ion_init[Z] = pyatomdb.atomdb.get_ionfrac(ionbalfile,
                      Z, te_init)
      ionfrac[Z][:,0] = ion_init[Z]

  # Deal with the ending index
  if (not pyatomdb.util.keyword_check(end_index)) or end_index < 0 \
     or end_index >= ncondi:
    end_index = ncondi-1

  # Some derived physical parameters
  temp_arr = conditions['kT']
  dens_arr = conditions['dens']
  radi_arr = conditions['R']
  velo_arr = conditions['velo']
  delr_arr     = np.zeros(ncondi, dtype=float)
  delr_arr[1:] = radi_arr[1:]-radi_arr[0:(ncondi-1)]
  meanv_arr     = np.zeros(ncondi, dtype=float)
  meanv_arr[0]  = velo_arr[0]
  meanv_arr[1:] = (velo_arr[1:]+velo_arr[0:(ncondi-1)])/2
  time_arr = delr_arr / meanv_arr
  meandens     = np.zeros(ncondi, dtype=float)
  meandens[0]  = dens_arr[0]
  meandens[1:] = (dens_arr[1:]+dens_arr[0:(ncondi-1)])/2
  tau_arr  = time_arr * meandens
  
  # Find the extreme values in the temperature array
  res1, _  = find_peaks(temp_arr)
  res2, _  = find_peaks(-temp_arr)
  extr_ind = np.sort(np.concatenate(([begin_index,end_index],res1,res2)))
  nextreme = len(extr_ind)
  
  for i in range(0,nextreme-1):
    iind = extr_ind[i]
    find = extr_ind[i+1]
        
    while iind < extr_ind[i+1]:
      taut = np.sum(tau_arr[(iind+1):(find+1)])
    
      # Derive the ionic fraction using the eletron temperature
      # at the beginning and the ending positions
      ionbal_l = {}
      ionbal_r = {}
      for Z in Zlist:
        ionbal_l[Z] = pyatomdb.apec.solve_ionbal_eigen(Z, temp_arr[iind], \
                      init_pop=ion_init[Z], tau=taut, teunit='keV')
      for Z in Zlist:
        ionbal_r[Z] = pyatomdb.apec.solve_ionbal_eigen(Z, temp_arr[find], \
                      init_pop=ion_init[Z], tau=taut, teunit='keV')
      # Compare the ionic fractions derived in one and two steps
      maxdiff = 0.0
      for Z in Zlist:
        # for iZ in range(0,Z+1):
        #   maxdiff = max([np.abs(ionbal_l[Z][iZ]-ionbal_r[Z][iZ]) / \
        #                    max([ionbal_l[Z][iZ],ionbal_r[Z][iZ],1e-5]), \
        #                  maxdiff])
          # print(maxdiff)
        maxdiff = max([max(np.abs(ionbal_l[Z]-ionbal_r[Z])), maxdiff])
    
      print(iind,find)
      print(ionbal_l[Z])
      print(ionbal_r[Z])
      print(maxdiff,taut/1e7)
      if maxdiff > diff_thres:
        if find == iind+1:
          if iind == begin_index:
            for Z in Zlist:
              ionfrac[Z] = ionbal_r[Z]
            ion_init = ionbal_r
            iind, find = find, min([find*2-iind+2, extr_ind[i+1]])
          else:
            print("The condition file at %d is too coarse!" % find)
            return -1
        else:
          find = int((iind+find)/2)
          continue
      else:
        accum_tau = np.cumsum(tau_arr[iind:(find+1)]) - tau_arr[iind]
        for Z in Zlist:
          for iZ in range(0,Z+1):
            ionfrac[Z][iZ,iind:(find+1)] = \
              np.interp(accum_tau,[0,taut], \
                   [ion_init[Z][iZ],ionbal_r[Z][iZ]])
        ion_init = ionbal_r
        iind, find = find, min([find*2-iind+2, extr_ind[i+1]])
        
  # # ionic fraction calculations
  # iind, find = begin_index, end_index
  # mind = int((iind+find) / 2)
  #
  # m_te = temp_arr[mind]
  # taut = np.sum(tau_arr)
  # ionbal_t = {}
  # for Z in Zlist:
  #   ionbal_t[Z] = pyatomdb.apec.solve_ionbal_eigen(Z, m_te, \
  #                 init_pop=ion_init[Z], tau=taut, teunit='keV')
  #
  # while iind < end_index:
  #   # Calculate the ionic fraction in one and two steps
  #   if mind > iind:
  #     # Tabled parameter
  #     i_te, m_te, f_te = temp_arr[[iind, mind, find]]
  #     # m_te = temp_arr[mind]
  #     # Derived parameter
  #     tau1 = np.sum(tau_arr[(iind+1):(mind+1)])
  #     # taut = tau1 + tau2
  #
  #     ionbal_1 = {}
  #     for Z in Zlist:
  #       ionbal_1[Z] = pyatomdb.apec.solve_ionbal_eigen(Z, m_te, \
  #                     init_pop=ion_init[Z], tau=tau1, teunit='keV')
  #
  #     if np.abs(i_te-m_te)/m_te > 1e-2 or np.abs(f_te-m_te)/m_te > 1e-2:
  #       find = mind
  #       mind = int((iind+find)/2)
  #       ionbal_t = ionbal_1
  #       continue
  #
  #     tau2 = np.sum(tau_arr[(mind+1):(find+1)])
  #     ionbal_2 = {}
  #     for Z in Zlist:
  #       ionbal_2[Z] = pyatomdb.apec.solve_ionbal_eigen(Z, m_te, \
  #                     init_pop=ionbal_1[Z], tau=tau2, teunit='keV')
  #     # Compare the ionic fractions derived in one and two steps
  #     maxdiff = 0.0
  #     for Z in Zlist:
  #       # for iZ in range(0,Z+1):
  #       #   maxdiff = max([np.abs(ionbal_t[Z][iZ]-ionbal_2[Z][iZ]) / \
  #       #                    max([ionbal_t[Z][iZ],1e-14]), maxdiff])
  #       maxdiff = max([max(np.abs(ionbal_t[Z]-ionbal_2[Z])),maxdiff])
  #     if maxdiff > diff_thres:
  #       find = mind
  #       mind = int((iind+find)/2)
  #       ionbal_t = ionbal_1
  #       continue
  #
  #   accum_tau = np.cumsum(tau_arr[iind:(find+1)]) - tau_arr[iind]
  #   for Z in Zlist:
  #     for iZ in range(0,Z+1):
  #       ionfrac[Z][iZ,iind:(find+1)] = \
  #         np.interp(accum_tau,[0,taut], \
  #              [ion_init[Z][iZ],ionbal_t[Z][iZ]])
  #   ion_init = ionbal_t
  #   iind, find = find, min([find*2-iind+2, end_index])
  #   mind = int((iind+find)/2)
    
  # Save calculated ionic fraction as pickle file
  tmp = open(outfilename,'wb')
  pickle.dump(ionfrac,tmp)
  tmp.close()
  return 0


#--------------------------------------------------------------------
#--------------------------------------------------------------------
#--------------------------------------------------------------------
def calc_nei_spectrum(Zlist, outfilename=False, condifile=False, \
  condi_index=False, ebins=False, ionfracfile=False, \
  linefile="$ATOMDB/apec_nei_line.fits", \
  cocofile="$ATOMDB/apec_nei_comp.fits", \
  dolines=True, appfile=False, rootpath=False):

  """
  calculate the NEI X-ray spectrum, based on the ionic fraction derived
  by the "calc_ionfrac" routine. The only input parameter is the index
  (array) of the elements.

  Parameters
  ----------
  Zlist: [int]
    list of element nuclear charges

  Keywords
  --------
  condifile: string or dictionary
    the ASCII file containing the physical condition array. can also
    pass in a structure read from the ASCII file;
  condi_index: [int]
    index array of at which condition position to derive the spectrum;
  ebins: [float]
    array of energy bins;
  ionfracfile: str or dictionary
    the file containing the ionic fraction at the all condition
    positions. Can also pass in the opened file, i.e.,
    "ionfrac = pickle.load(open('tionfrac.pkl','rb'))
  linefile: str or HDUList
    The file containing all the line emission. Defaults to
    "$ATOMDB/apec_nei_line.fits". Can also pass in the opened file,
    i.e. "linefile = pyatomdb.pyfits.open('apec_nei_line.fits')";
  cocofile: str or HDUList
    The file containing all the continuum emission. Defaults to
    "$ATOMDB/apec_nei_comp.fits". Can also pass in the opened file,
    i.e. "cocofile = pyatomdb.pyfits.open('apec_nei_comp.fits')";
  outfilename: str
    the name of the output pickle file;
  appfile: str
    If set, the pickle file will be update, otherwise a new pickle
    named "tspec_Element.String.pkl" will be created.

  Returns
  -------
  No return, but the pickle file is created/updated with derived NEI
  spectra at condition positions specified by condi_index.

  """

  # System parameter
  if not pyatomdb.util.keyword_check(rootpath):
    rootpath = os.getcwd()+'/'

  # Output file name
  if not pyatomdb.util.keyword_check(outfilename):
    if pyatomdb.util.keyword_check(appfile):
      outfilename = appfile
    else:
      outfilename = 'tspec_'
      for Z in Zlist:
        outfilename += pyatomdb.atomic.Ztoelsymb(Z)
      outfilename += '.pkl'

  # Deal with the element list
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
      conditions = ascii.read(confile)
    elif isinstance(condifile, astropy.table.table.Table):
      conditions = condifile
    else:
      print("Unknown data type for condition file. Please pass a " \
        "string or an ASCIIList")
      return -1
  ncondi = len(conditions)
  if not pyatomdb.util.keyword_check(condi_index):
    condi_index = range(0,ncondi)
  else:
    if max(condi_index) >= ncondi:
      print("*** Supplied condition index is illegal. Exiting. ***")
      return -1

  # Set up spectral bins
  if not pyatomdb.util.keyword_check(ebins):
    mineng = 0.1
    maxeng = 2.0
    nbins = 190
    ebins = np.linspace(mineng,maxeng,nbins)
  else:
    nbins = len(ebins)

  # Read the ionic fraction
  if not pyatomdb.util.keyword_check(ionfracfile):
    ionfracfile = 'tionfrac_'
    for Z in Zlist:
      ionfracfile += pyatomdb.atomic.Ztoelsymb(Z)
    ionfracfile += '.pkl'
  if isinstance(ionfracfile, str):
    ionfile = os.path.expandvars(rootpath+ionfracfile)
    if not os.path.isfile(ionfile):
      print("*** ERROR: no such ionic fraction file %s. Exiting ***" \
        % (ionfile))
      return -1
    ionfrac = pickle.load(open(ionfile,'rb'))
  elif isinstance(ionfracfile, dict):
    ionfrac = ionfracfile
  else:
    print("Unknown data type for ionic fraction file. Please pass a" \
      "string or an ASCIIList.")
    return -1

  #Read the emissivity files, if "linefile" and "cocofile" are strings
  if isinstance(linefile, str):
    lfile = os.path.expandvars(linefile)
    if not os.path.isfile(lfile):
      print("*** ERROR: no such file %s. Exiting ***" %(lfile))
      return -1
    ldat = pyfits.open(lfile)
  elif isinstance(linefile, pyfits.hdu.hdulist.HDUList):
    # no need to do anything, file is already open
    ldat = linefile
  else:
    print("Unknown data type for linefile. Please pass a string or an HDUList")
    return -1
  if isinstance(cocofile, str):
    cfile = os.path.expandvars(cocofile)
    if not os.path.isfile(cfile):
      print("*** ERROR: no such file %s. Exiting ***" %(cfile))
      return -1
    cdat = pyfits.open(cfile)
  elif isinstance(cocofile, pyfits.hdu.hdulist.HDUList):
    # no need to do anything, file is already open
    cdat = cocofile
  else:
    print("Unknown data type for cocofile. Please pass a string or an HDUList")
    return -1
  
  # Find HDU with kT closest to desired kT in given line or coco file
  # ite = pyatomdb.spectrum.get_index(trans_te)
  # Get HDU with kT closest to-but-less than desired kT
  trans_te = conditions['kT']
  fl_ind = np.log10(trans_te/pyatomdb.const.KBOLTZ)*10-40+2
  ite    = np.int16(fl_ind)
  fact1  = 1.0-(fl_ind-ite)
  fact2  = fl_ind-ite

  # The final result - spec_total
  spec_total = {}
  spec_total['ebins'] = ebins
  for Z in Zlist:
    spec_total[Z] = np.zeros([ncondi,nbins],dtype=float)

  # If the keyword "appfile" is set, append it.
  if pyatomdb.util.keyword_check(appfile):
    old_spec = pickle.load(open(rootpath+appfile,'rb'))
    # mod_spec = {}
    for Z in Zlist:
      # mod_spec[Z] = np.zeros([ncondi,nbins],dtype=float)
      ncondi_old_spec = len(old_spec[Z][:,0])
      if ncondi_old_spec < ncondi:
        spec_total[Z][0:ncondi_old_spec,:] = old_spec[Z]
      else:
        spec_total[Z] = old_spec[Z][0:ncondi,:]

  for l in condi_index:
    print('For Zone-%03d:...' %l, end='', flush=True)
    for Z in Zlist:
      one_spec = np.zeros(nbins,dtype=float)
      for iZ in range(1,Z+2):
        spec1 = pyatomdb.spectrum.make_ion_spectrum(ebins, ite[l], Z, \
                  iZ, nei=True, dummyfirst=True, linefile=ldat, \
                  cocofile=cdat, dolines=dolines)
        spec2 = pyatomdb.spectrum.make_ion_spectrum(ebins, ite[l]+1, Z,
                  iZ, nei=True, dummyfirst=True, linefile=ldat, \
                  cocofile=cdat, dolines=dolines)
        tspec = spec1*fact1[l]+spec2*fact2[l]
        one_spec += ionfrac[Z][iZ-1,l]*tspec

      spec_total[Z][l,:] = one_spec
    print("  finished.")

  # Save derived spectrum as pickle file
  tmp = open(rootpath+outfilename,'wb')
  pickle.dump(spec_total,tmp)
  tmp.close()
