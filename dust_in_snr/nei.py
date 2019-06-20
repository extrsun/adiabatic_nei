# -*- coding: utf-8 -*-
"""
Created on Wed Jul 12 20:53:21 2017

@author: zehaojin
"""

from scipy.integrate import odeint
from ctypes import *
import matplotlib.pyplot as plt
from numpy import *
import pyatomdb



def calc_elem_ion_rec_rate(Z,Te,TkeV):
#调用pyatomdb获取每个元素的电离率和复合率
    print('element=%d Temperature=%f keV\n'% (Z, TkeV))
    ###create a list called rate
    rate=[[0.0]*(Z+1),[0.0]*(Z+1)]#定义数据结构rate=[27个[0.0]，27个[0.0]]
    #Z+1表示Fez，因为o是尘埃里的fe原子，fe1表示气体的原子；所以fe27表示的是电子全部剥离的26价铁离子.
    #Z具体的赋值是在ionfrac_data.py里面调用的，含义就是原子序数，这个是不会变的，就是后面的这个26，yy=nei.calc_nei(26, taulist, fe, Te_K,telist[ite])
    for i in range(1,Z+1):#就是1到27，中性铁原子到26价的铁离子
        ionrec_rate=pyatomdb.atomdb.get_ionrec_rate(Te, 'irdat_in', lvdat_in=False, Te_unit='K', lvdatp1_in=False, ionpot=False, separate=False, Z=Z, z1=i, settings=False, datacache=False, extrap=True)
        ###Te：electron temperature in K (default), eV, or keV
        ###returns:(ionization rate coeff., recombination rate coeff.) in cm^3 s^-1
        rate[0][i]=ionrec_rate[0] #ionization rate coeff
        rate[1][i]=ionrec_rate[1] #recombination rate coeff
    print('ionization rate_under standard nei model(0-27 ion):\n',rate[0],'\n')
    print('recombination rate_under standard nei model(0-27 ion):\n',rate[1],'\n')
    return  rate 
    #所以最后的这个数组就是电离率1到27，中性铁原子到26价的铁离子，以及复合率1到27，中性铁原子到26价的铁离子



def calc_nei(Z,timegrid,init,TeK,TkeV):
	dustZele=init[0]
	#initial ion_fractions ,fe=numpy.zeros(28,dtype=float);fe[0]=dust_init;fe[1]=1.0-dust_init
	#so, dustZele=init[0]=fe[0]=dust_init=input=95/100=0.95
	rate=calc_elem_ion_rec_rate(Z,TeK,TkeV)
	def func(Zele,tau):#定义微分方程组
		y=zeros(Z+2,dtype=float)#定义函数y有28个值的数组。
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
TK=1e7
tau=linspace(1,1e13,1000001)
ele=zeros(28,dtype=float)
ele[0]=0
ele[1]=1.0
y=calc_nei(26,tau,ele,TK)
fig1=plt.figure('fig1')
for z in range(0,28):
	a=[]
	for k in y:
		a.append(k[z])
	plt.plot(tau,a)
plt.loglog()
plt.xlim(1e9,1e13)
plt.ylim(1e-6,1)
plt.show()
'''
