import nei,pylab, pyatomdb,numpy,pickle
from scipy import *

K2keV=8.62e-8
#能量单位转换，k to kev;1 kelvin [K] = 8.61732814974056E-08 kiloelectron-volt [keV]

telist=numpy.array([1.,2.,4.,10.])
#设定的尘埃的环境温度，分别是1，2，4，10kev，可根据需要修改

E_keV=[0,6.4077,6.3915]
#For isolated Fe atoms, the ﬂuorescent Fe Ka line consists of two components,K α 1 =6.404keV 和 K α 2 =6.391keV

yie=[0,0.2026,0.1013] 
#florescent yields;K-shell ionization rates must be multiplied by the ﬂuorescent yields for production of Ka photons (equal to 0.30 for Fe)


#################################################
version='3.0.9'
#################################################