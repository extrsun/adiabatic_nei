#calculate ionic fraction
#input initial dust percentage dust_perc
#initial condition: dust x, atomic Fe 1-x, other ions 0.
#store results (dat) in ionfrac_dust_perc_v_version.pkl
#dat is a dictionary,
#      dat['telist']:  electron temerature
#			 dat['taulist']: tau g rids, from 1 to 1e13, linearly distributed
#			 dat['ionfrac']: three-dimension list, ionfrac[Te,tau,ion], ion from 0 (dust) to 27 (FeXXVII or Fe26+)

from constants import *
from datetime import datetime

taulist = numpy.linspace(1,1e13,100001)
ionfrac = numpy.zeros([len(telist),len(taulist),28],dtype=float)

dust_perc=input('please enter the percentage of dust, e.g. 95:\n') 
print('percentage of dust=', dust_perc)
print('working on it...')
dust_init=float(dust_perc)/100.

now1 = datetime.now().hour*3600. + datetime.now().minute*60. + \
       datetime.now().second + datetime.now().microsecond/1e6
for ite in range(len(telist)):
	fe=numpy.zeros(28,dtype=float)
	fe[0]=dust_init
	fe[1]=1.0-dust_init
	Te_K=telist[ite]/K2keV 
	#convert kev to kelven temprature.	
	yy=nei.calc_nei(26, taulist, fe, Te_K,telist[ite])
	'''this is only for reference. actually it is copied from nei.py.
	def calc_nei(Z,timegrid,init,TeK,TkeV):
	dustZele=init[0]
	rate=calc_elem_ion_rec_rate(Z,TeK,TkeV)
	def func(Zele,tau):
		y=zeros(Z+2,dtype=float)
		tau_13=(1e-13)*tau
		f_tau=1.0-3.2*tau_13**0.5+3.0*tau_13-tau_13**2+0.2*tau_13**3
		df_tau=-1.6/(tau_13**0.5)+3.0-2*tau_13+0.6*tau_13**2
		
		y[0]=dustZele*(df_tau*10**(-13))
		y[1]=dustZele*(-df_tau*10**(-13))-rate[0][1]*Zele[1]+rate[1][1]*Zele[2]
		for k in range(2,Z+1):
			y[k]=-(rate[0][k]+rate[1][k-1])*Zele[k]+rate[0][k-1]*Zele[k-1]+rate[1][k]*Zele[k+1]
		y[Z+1]=-rate[1][Z]*Zele[Z+1]+rate[0][Z]*Zele[Z]
		
		return y
	result=odeint(func,init,timegrid)
	return result

	'''
	for itau in range(len(taulist)):
		ionfrac[ite,itau,:]=yy[itau]
print('dust inclued ionfrac list[Te,tau,ion]=\n',ionfrac)

now2 = datetime.now().hour*3600. + datetime.now().minute*60. + \
       datetime.now().second + datetime.now().microsecond/1e6
print("Time Consuming:%7.2f sec." % (now2-now1))

dat={}
dat['telist']=telist
dat['taulist']=taulist
dat['ionfrac']=ionfrac
pickle.dump(dat,open('data_package/ionfrac_%s_v%s.pkl' %(str(dust_perc),version),'wb'))
###pickle.dump(dat,open('data_package/ionfrac_'+str(dust_perc)+'.pkl','wb'))
