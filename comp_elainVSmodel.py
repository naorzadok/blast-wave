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
W15_l  = np.array([2.1777312115225382, 4.300583148173627, 6.788872298048224, 8.214338980815457, 9.78481713569536, 13.366887021525171, 17.487472255192074, 22.08963925769262, 27.047272365150658])
W16_l = np.array([2.0099448159652655, 3.8616304616834367, 6.0017247326066645, 7.2268213435406, 8.577999547773649, 11.664230325132511, 15.216815094569808, 19.182325129368397, 23.447555070072138])
W17_l  = np.array([1.8730000384218584, 3.519578964454646, 5.403686489984404, 6.481947664699828, 7.671823706440207, 10.389957301074253, 13.515075278859827, 16.995322765673777, 20.72687559490444])
W18_l  = np.array([1.7589620627028084, 3.2441126362937047, 4.929205921488902, 5.892669571433321, 6.955487129659217, 9.37999040501026, 12.158743567029541, 15.240494083364133, 18.528655525151002])
W19_l  = np.array([1.6608343647548551, 3.0130689074920487, 4.534234786704071, 5.402025040288163, 6.357924072799259, 8.53187122726117, 11.010374364833378, 13.742251384568663, 16.636890948488958])
W20_l  = np.array([1.5758760243319276, 2.815652964021892, 4.1964967701014455, 4.981097369401946, 5.8430466509480565, 7.793742614497345, 10.000741586596885, 12.412653092906092, 14.944116159529102])
E_l  = np.array([2,4,6,7,8,10,12,14,16])
Edot = np.array([17.274618317056497,23.728294177732373,34.111531502363086,28.013687149637473,27.546570548130653,28.01766980360358,22.917351334547227,17.523004124482284,14.299547152187518])
E2_l  = np.array([2,4,6,8,10,12,14,16,18])
x_ind = 0.175
x15_ax  = W15_l/x_ind
x16_ax  = W16_l/x_ind
x17_ax  = W17_l/x_ind
x18_ax  = W18_l/x_ind
x19_ax  = W19_l/x_ind
x20_ax  = W20_l/x_ind

eff  = np.linspace(2,20,19)
res  = f_el(eff)


plt.plot(E_l,x15_ax,label = 'overdrive 1.5')
plt.plot(E_l,x16_ax,label = 'overdrive 1.6')
plt.plot(E_l,x17_ax,label = 'overdrive 1.7')
plt.plot(E_l,x18_ax,label = 'overdrive 1.8')
plt.plot(E_l,x19_ax,label = 'overdrive 1.9')
plt.plot(E_l,x20_ax,label = 'overdrive 2.0')
plt.plot(eff,res,label = 'Lu et al')
plt.plot(E2_l,Edot,'o',label = 'Lu et al')
plt.grid()
plt.ylabel('$\u03BB/{x_d}$',size=18)
plt.xlabel('effective activation energy',size=18)
plt.legend(loc='upper right',fontsize = 16)
plt.xlim(0,20)
plt.ylim(0,70)
plt.xticks(size=16)
plt.yticks(size=16)
plt.show()