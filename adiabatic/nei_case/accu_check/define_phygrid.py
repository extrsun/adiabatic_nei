"""
Purpose: pre-define three physical grids for accuracy comparison.

Written by sw, Jun 19, 2019, initial release
"""

import pickle, os, math
import numpy as np
import scipy.constants as cons

def def_phygrid(nbins, outname="phys_grid.pkl", log=False, \
                   rootpath=os.getcwd()+'/'):
  """
  Sub-routine that defines the physical grid.
  
  Parameters
  ----------
  nbins: int
    number of the physical grids
  
  Keywords
  --------
  outname: str
    name of the output .pkl file (def: "phys_grid.pkl")
  log: boolean
    whether draw the grid in logarithmic scale or not
  rootpath: str
    root path that places the output pickle file (def: current path)
  
  Return
  ------
  No return, but the pickle file with outname is created.
  
  """
  
  # Physical parameters defined in Ji+2006, Sec.3.4
  parsec = 3.0856780e+18
  mu  = 1.4
  u0  = 1e3 #in km/s
  T0  = 5e6 #in K
  cs0 = np.sqrt( 5*cons.k*1e7*T0 / (3*mu*cons.m_p*1e3) ) / 1e5 #in km/s
  Tc  = T0 + mu*cons.m_p*1e3*(u0*1e5)**2 / (5*cons.k*1e7) #in K
  r0  = 0.3 #in pc
  mdot0 = 3e-5 * 1.989e+33 / cons.year #in g/s
  c1  = mdot0 / (4*math.pi) #constant coefficient for \rho*r**2*u
  nH0 = c1 / ( u0*1e5 * (r0*parsec)**2 * mu*cons.m_p*1e3 ) #in cm-3
  ne0 = nH0 * (mu+1)/2
  
  if log:
    Tval = 10**( -np.linspace(0,nbins,nbins+1) / (nbins/2.5) )*T0
  else:
    Tval = np.linspace(nbins,1,nbins) / nbins * T0
  nTval  = len(Tval)
  
  # Derive the radius array and the other physical parameters
  Tr_val = r0 / ( (Tval/T0)**3 * (Tc-Tval)/(Tc-T0) )**0.25 #in pc
  u_app  = np.sqrt( 3*cs0**2 + u0**2 - 3*cs0**2/(Tr_val/r0)**(4/3.) )
  ne_app = u0/u_app * (r0/Tr_val)**2 * ne0
  
  # Save the physical grid as pickle file
  outdat = {}
  outdat['R'] = Tr_val*parsec
  outdat['velo'] = u_app*1e5
  outdat['dens'] = ne_app
  outdat['kT']   = Tval*cons.k*1e7 / (cons.eV*1e3*1e7)
  tmp = open(rootpath+outname,'wb')
  pickle.dump(outdat,tmp)
  tmp.close()
  
  
if __name__ == "__main__":
  rootpath = os.getcwd()+'/'  
  def_phygrid(300, outname='adia_exp.phy_300_log.pkl', log=True, \
    rootpath=rootpath)
  def_phygrid(3000, outname='adia_exp.phy_3000_log.pkl', log=True, \
    rootpath=rootpath)
  def_phygrid(300, outname='adia_exp.phy_300_linear.pkl', \
    rootpath=rootpath)