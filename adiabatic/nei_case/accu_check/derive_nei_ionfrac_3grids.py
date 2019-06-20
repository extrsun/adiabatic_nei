import nei_physcalc
import pickle, os
import pyatomdb
from datetime import datetime
import multiprocessing as mp

#system parameters
rootpath = os.getcwd()+'/'
# Zlist = [1,2,6,7,8,10,12,14,16,18,20,26,28]
# Zlist = [6,7,8,10,12,14]
# Zlist = [1,2,16,18,20,26,28]
Zlist = [6]
modname = ['300_log', '300_linear', '3000_log']
nmod = 3

# calculate the ionic fraction in parallel way-------------
now1 = datetime.now().hour*3600. + datetime.now().minute*60. + \
       datetime.now().second + datetime.now().microsecond/1e6
pool = mp.Pool(mp.cpu_count())
res  = pool.starmap(nei_physcalc.calc_nei_ionfrac, \
         [(Zlist, 'adia_exp.phy_'+modname[Z]+'.pkl', \
         'tionfrac_C_'+modname[Z]) for Z in range(0,nmod)])
now2 = datetime.now().hour*3600. + datetime.now().minute*60. + \
       datetime.now().second + datetime.now().microsecond/1e6
print("Time Consuming:%7.2f sec." % (now2-now1))

# for Z in range(0,nmod):
#   now1 = datetime.now().hour*3600. + datetime.now().minute*60. + \
#          datetime.now().second + datetime.now().microsecond/1e6
#   nei_physcalc.calc_nei_ionfrac(Zlist, condifile='adia_exp.phy_'+ \
#     modname[Z]+'.pkl', outfilename='tionfrac_C_'+modname[Z]+'.pkl')
#   now2 = datetime.now().hour*3600. + datetime.now().minute*60. + \
#          datetime.now().second + datetime.now().microsecond/1e6
#   print("# Time Consuming:%7.2f sec for model %s" % (now2-now1, modname[Z]))
