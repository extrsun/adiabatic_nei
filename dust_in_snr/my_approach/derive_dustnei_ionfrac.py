#Purpose: derive the NEI ionic fraction in the case of dusty SNR

# Import used modules
import dustnei_physcalc
import pickle, os
import numpy as np
from datetime import datetime
import multiprocessing as mp

def derive_dustnei_ionfrac_module(log=False):
  # System parameters
  rootpath = os.getcwd() + '/'
  if log:
    rootpath += 'log/'
  else:
    rootpath += 'linear/'
  Zlist = [26]

  # Define the condition dictionary, based on Zehao's setting
  kTlist  = np.array([1.,2.,4.,10.]) # in keV
  if log:
    taulist = np.zeros(1201, dtype=float) + 1e10 # in cm-3 s
    taulist[0:301] = 10**(np.linspace(9,10,301))
  else:
    taulist = np.zeros(1001, dtype=float) + 1e10 # in cm-3 s
  # taulist = 10**(np.linspace(9,11,1001)) # in cm-3 s
  # taulist = np.array([10**(np.linspace(9,10,301)), \
  #                    np.zeros(900) + 1e10]) # in cm-3 s
  accum_tau = np.cumsum(taulist) - taulist[0]
  n_kT  = len(kTlist)
  n_tau = len(taulist)
  nelist = np.zeros(n_tau, dtype=float) + 1.0  # in cm-3
  v_list = np.zeros(n_tau, dtype=float) + 1e10 # in cm/s
  r_list = accum_tau / nelist * v_list # in cm

  conditions = {}
  conditions['R'] = r_list
  conditions['velo'] = v_list
  conditions['dens'] = nelist
  conditions['kT'] = np.zeros(n_tau, dtype=float) + kTlist[0]

  # Derive the change of the neutral iron
  dustZele  = 0.95
  tau_13 = accum_tau/1e13
  f_tau=1.0-3.2*tau_13**0.5+3.0*tau_13-tau_13**2+0.2*tau_13**3
  init_ionfrac = np.zeros([27, n_tau], dtype=float)
  init_ionfrac[0,:] = (1-dustZele) + dustZele * (1-f_tau)
  init_file = {}
  for Z in Zlist:
    init_file[Z] = init_ionfrac

  now1 = datetime.now().hour*3600. + datetime.now().minute*60. + \
         datetime.now().second + datetime.now().microsecond/1e6
  pool = mp.Pool(mp.cpu_count())
  res  = pool.starmap(dustnei_physcalc.calc_nei_ionfrac, \
           [(Zlist, conditions, 'tionfrac_Fe_%04.1fkeV.pkl'%kTlist[M], \
            init_file, kTlist[M]) for M in range(0,n_kT)])
  # for M in range(0,n_kT):
  #   condition['kT'] = np.zeros(n_tau, dtype=float) + kTlist[M] # in keV
  #   dustnei_physcalc.calc_nei_ionfrac(Zlist, condifile=condition, \
  #     outfilename='tionfrac_Fe_%fkeV.pkl'%kTlist[M], \
  #     init_file=init_file, rootpath=rootpath)
  now2 = datetime.now().hour*3600. + datetime.now().minute*60. + \
         datetime.now().second + datetime.now().microsecond/1e6
  print("# Time Consuming:%7.2f sec." % (now2-now1))

  # Combine derived ionic fraction
  ionfrac = np.zeros([n_kT, 28, n_tau], dtype=float)
  for M in range(0,n_kT):
    res = pickle.load(open(rootpath+'tionfrac_Fe_%04.1fkeV.pkl'%kTlist[M], 'rb'))
    ionfrac[M,0,:] = dustZele * f_tau
    for iZ in range(0,Z+1):
      ionfrac[M,iZ+1,:] = res[Z][iZ,:] * init_ionfrac[0,:]

  dat={}
  dat['telist']  = kTlist
  dat['taulist'] = accum_tau
  dat['ionfrac'] = ionfrac
  tmp = open(rootpath+'dustsnr_ionfrac_Fe.pkl', 'wb')
  pickle.dump(dat, tmp)
  tmp.close()

if __name__ == '__main__':
  derive_dustnei_ionfrac_module(log=True)
  derive_dustnei_ionfrac_module(log=False)