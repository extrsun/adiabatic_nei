"""
The physcalc module contains the routine that gathering the ionic
fraction and the X-ray spectrum of CIE based on the physical
conditions in an ASCII file.

V0.1 - Wei Sun, May 28, 2019: initial release
V0.2 - Wei Sun, Jun 01, 2019: standardize for github uploading
V0.3 - Wei Sun, Jun 03, 2019: add keyword "dolines" to calculate
       pure continuum spectrum
"""

import pickle, os
import numpy as np
import pyatomdb
from astropy.io import ascii
import astropy

def deri_cie_ionfrac(Zlist, condifile='adia.exp_phy.info', \
  condi_index=False, appfile=False, outfilename=False, rootpath=False):

  """
  Derive the CIE ionic fraction based on the physical conditions in
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
  condi_index: [int]
    index array of at which condition position to derive the spectrum;
  appfile: str or dictionary of ionic fraction
    the pickle file that the new calculation will be appended into.
    Could be the dictionary of loaded from the pickle file;
  outfilename: str
    the name of output pickle file recording the ionic fraction.
    The name of the output file is adopted as following sequence:
      1. specified by <outfilename>;
      2. adopted from <appfile>, if applicable;
      3. "tionfrac_Element.List.pkl".

  Returns
  -------
  No return, but a pickle file containing derived CIE ionic fraction at
  specified condition positions is created/updated.
  """

  # System parameter
  atomdbpath = os.environ['ATOMDB']
  ionbalfile = atomdbpath+'APED/ionbal/v3.0.7_ionbal.fits'
  if not pyatomdb.util.keyword_check(rootpath):
    rootpath = os.getcwd()+'/'

  NZlist = len(Zlist)

  # Output file name
  if not pyatomdb.util.keyword_check(outfilename):
    if pyatomdb.util.keyword_check(appfile):
      if isinstance(appfile, str):
        outfilename = appfile
    else:
      outfilename = 'tciefrac_'
      for Z in Zlist:
        outfilename += pyatomdb.atomic.Ztoelsymb(Z)
      outfilename += '.pkl'

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
      return -1

  te_arr = conditions['kT']/pyatomdb.const.KBOLTZ #in K

  ionfrac = {}
  for Z in Zlist:
    ionfrac[Z] = np.zeros([Z+1,ncondi],dtype=float)

  for l in condi_index:
    print('For Zone-%03d...' % l, end='')
    for Z in Zlist:
      ionfrac[Z][:,l] = pyatomdb.atomdb.get_ionfrac(ionbalfile, \
                          Z, te_arr[l])
    print('finished.')

  # Save calculated ionic fraction as pickle file
  tmp = open(outfilename,'wb')
  pickle.dump(ionfrac,tmp)
  tmp.close()
  return 0

def deri_cie_spectrum(Zlist, condifile="adia.exp_phy.info", \
  condi_index=False, ebins=False, outfilename=False, appfile=False, \
  linefile="$ATOMDB/apec_line.fits",cocofile="$ATOMDB/apec_coco.fits", \
  dolines=True, rootpath=False):

  """
  calculate the CIE X-ray spectrum, based on the ionic fraction derived
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
    "ionfrac = pickle.load(open('tionfrac_cie.pkl','rb'))
  linefile: str or HDUList
    The file containing all the line emission. Defaults to
    "$ATOMDB/apec_nei_line.fits". Can also pass in the opened file,
    i.e. "linefile = pyatomdb.pyfits.open('apec_line.fits')";
  cocofile: str or HDUList
    The file containing all the continuum emission. Defaults to
    "$ATOMDB/apec_nei_comp.fits". Can also pass in the opened file,
    i.e. "cocofile = pyatomdb.pyfits.open('apec_coco.fits')";
  outfilename: str
    the name of the output pickle file;
  appfile: str
    If set, the pickle file will be update, otherwise a new pickle
    named "tspec_Element.String.pkl" will be created.

  Returns
  -------
  No return, but the pickle file is created/updated with derived CIE
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
      outfilename = 'tcspec_cie.pkl'

  # Deal with the element list
  NZlist = len(Zlist)

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

  # Set up spectral bins
  if not pyatomdb.util.keyword_check(ebins):
    mineng = 0.1
    maxeng = 2.0
    nbins = 190
    ebins = np.linspace(mineng,maxeng,nbins)
  else:
    nbins = len(ebins)

  # Deal with the appendix file
  spec_total = np.zeros([ncondi,nbins],dtype=float)
  if pyatomdb.util.keyword_check(appfile):
    old_spec = pickle.load(open(rootpath+appfile,'rb'))
    ncondi_old_spec = len(old_spec[:,0])
    if ncondi_old_spec < ncondi:
      mod_spec = np.zeros([ncondi,nbins],dtype=float)
      mod_spec[0:ncondi_old_spec,:] = old_spec
    else:
      mod_spec = old_spec
    sind = np.where(mod_spec[:,1] != 0)
    if len(sind[0]) > 0:
      for i in sind[0]:
        spec_total[i,:] = mod_spec[i,:]

  # Deal with the condition index
  if not pyatomdb.util.keyword_check(condi_index):
    if pyatomdb.util.keyword_check(appfile):
      sind = np.where(mod_spec[:,1] == 0)
      if len(sind[0]) > 0:
        condi_index = sind[0]
      else:
        print("*** Nothing to do! ***")
        return 0
    else:
      condi_index = range(0,ncondi)
  else:
    if max(condi_index) >= ncondi:
      print("*** Supplied condition index is illegal. Exiting. ***")
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
  
  # Tabled parameter
  trans_te = conditions['kT'] #temperature of zone
  trans_R  = conditions['R'] #radius of zone
  # # Find HDU with kT closest to desired kT in given line or coco file
  # ite = pyatomdb.spectrum.get_index(trans_te)
  # Get HDU with kT closest to-but-less than desired kT
  fl_ind = np.log10(trans_te/pyatomdb.const.KBOLTZ)*10-40+2
  ite    = np.int16(fl_ind)
  fact1  = 1.0-(fl_ind-ite)
  fact2  = fl_ind-ite

  # The radius/zone cycle
  for l in condi_index:
    print('For Reg-%03d: R=%10.3e, kT=%4.2f keV:...' % \
      (l, trans_R[l], trans_te[l]), end='')

    # Derive the spectrum for one zone
    spec1 = pyatomdb.spectrum.make_spectrum(ebins, ite[l],
              elements = Zlist, dummyfirst=True, \
              linefile=ldat, cocofile=cdat, \
              dolines=dolines)
    spec2 = pyatomdb.spectrum.make_spectrum(ebins, ite[l]+1,
              elements = Zlist, dummyfirst=True, \
              linefile=ldat, cocofile=cdat, \
              dolines=dolines)
    tspec = spec1*fact1[l]+spec2*fact2[l]
    spec_total[l,:]=tspec
    print('  Finished.')

  # Save test spectra as pickle file
  tmp = open(rootpath+outfilename,'wb')
  pickle.dump(spec_total,tmp)
  tmp.close()
