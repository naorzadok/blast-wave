from re import X
from typing import List
import numpy as np
from numpy.core.defchararray import array
from numpy.core.fromnumeric import argmin, size
from scipy.integrate.quadpack import quad
import scipy.optimize as opt
import matplotlib.pyplot as plt
import math
import sympy as sym
from scipy.misc import derivative as dev
from scipy.optimize import zeros
from scipy.optimize.moduleTNC import minimize
from sympy.core import symbol
from sympy.utilities.lambdify import lambdify
from scipy.integrate import *

# units in centimeters
f_el = lambda x: -5.34e-4*x**4+3.63e-2*x**3-9.41e-1*x**2+8.96*x+2.
W15_l  = np.array([1.4220345115030826, 2.51017907715519, 3.566147580141776, 4.709214807444818])
W16_l = np.array([1.3142224876525281, 2.2531968317912887, 3.1461967177948926, 4.112009273030238])
W17_l  = np.array([1.2252184706039986, 2.049211062377554, 2.8184471718606416, 3.6474537596674668])
W18_l  = np.array([1.1482449206017877, 1.8778376399544938, 2.5450225949308867, 3.2582809851323336])
W19_l  = np.array([1.0796957551063804, 1.7278512850831267, 2.3054791931654424, 2.9144063919182126])
W20_l  = np.array([1.0152629412605747, 1.5891099772378587, 2.0834670755095757, 2.5934593032141677])
E_l  = np.array([2,4,6,8])
E_res = np.array([3.6,3.3,3.48774919556346])
E2_l  = np.array([4,6,8])
x_ind = 0.175
x15_ax  = W15_l/x_ind
x16_ax  = W16_l/x_ind
x17_ax  = W17_l/x_ind
x18_ax  = W18_l/x_ind
x19_ax  = W19_l/x_ind
x20_ax  = W20_l/x_ind
Edot    = E_res/x_ind

eff  = np.linspace(2,20,19)
res  = f_el(eff)


plt.plot(E_l,x15_ax,label = 'overdrive 1.5')
plt.plot(E_l,x16_ax,label = 'overdrive 1.6')
plt.plot(E_l,x17_ax,label = 'overdrive 1.7')
#plt.plot(E_l,x18_ax,label = 'overdrive 1.8')
#plt.plot(E_l,x19_ax,label = 'overdrive 1.9')
#plt.plot(E_l,x20_ax,label = 'overdrive 2.0')
#plt.plot(eff,res,label = 'Lu et al_not quite')
plt.plot(E2_l,Edot,'o',label = 'gamma =1.5')
plt.grid()
plt.ylabel('$\u03BB/{x_d}$',size=18)
plt.xlabel('effective activation energy',size=18)
plt.legend(loc='upper right',fontsize = 16)
plt.xlim(0,20)
plt.ylim(0,70)
plt.xticks(size=16)
plt.yticks(size=16)
plt.show()