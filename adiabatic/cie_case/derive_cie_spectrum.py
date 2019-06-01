#Purpose: derive the CIE spectrum

# Import used modules
import cie_physderi
import pyatomdb
import pickle, os
from datetime import datetime
from astropy.io import ascii
import numpy as np

#system parameters
rootpath = os.getcwd()+'/'
Zlist = [1,2,6,7,8,10,12,14,16,18,20,26,28]

#Read the condition file
confile = rootpath+'adia.exp_phy.info'
conditions = ascii.read(confile)
ncondi = len(conditions)
condi_index = range(0,ncondi)

# set up spectral bins
mineng = 0.1
maxeng = 2.0
nbins = 190
ebins = np.linspace(mineng,maxeng,nbins)

# pre-open the emissivity files
linefile  = \
  pyatomdb.pyfits.open(os.path.expandvars('$ATOMDB/apec_line.fits'))
cocofile  = \
  pyatomdb.pyfits.open(os.path.expandvars('$ATOMDB/apec_coco.fits'))

# spectrum calculation
now1 = datetime.now().hour*3600. + datetime.now().minute*60. + \
       datetime.now().second + datetime.now().microsecond/1e6
cie_physderi.deri_cie_spectrum(Zlist, condifile=conditions, \
               condi_index=condi_index, ebins=ebins, \
               linefile=linefile, cocofile=cocofile, \
               outfilename="tcspec_cie.pkl")
now2 = datetime.now().hour*3600. + datetime.now().minute*60. + \
       datetime.now().second + datetime.now().microsecond/1e6
print("Time Consuming:%7.2f sec." % (now2-now1))
