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
W15_l  = np.array([2.014601666382274, 3.9027358298764057, 6.037769732226555, 8.574315656854441, 11.549781216008206, 15.0055686376287, 18.84094884652358, 23.00915647360085])
W16_l = np.array([1.8605624801486362, 3.5060203402249948, 5.339494521048837, 7.518097151951526, 10.07822824495078, 13.053484123135435, 16.352140659211138, 19.930127963745647])
W17_l  = np.array([1.7344630623786976, 3.195830623920909, 4.806609313013051, 6.72047398800872, 8.969274750305896, 11.578720068758706, 14.463038501148086, 17.57993786366952])
W18_l  = np.array([1.629966447585842, 2.9464825731578608, 4.383652201584381, 6.088505926873556, 8.086918333646977, 10.39715040707943, 12.937638397204305, 15.667285090269019])
W19_l  = np.array([1.5403155217868993, 2.7372906294697352, 4.030608946279873, 5.55882761495993, 7.341293787903419, 9.38926725761701, 73.5537983625194, 14.00658150529158])
W20_l  = np.array([1.4608004967561503, 2.5546250022193995, 3.7220975864192245, 5.092261470576819, 6.677810509318696, 8.483451268303698, 10.433695101069407, 12.489316429932776])
E_l  = np.array([2,4,6,8,10,12,14,16])
E_res = np.array([3.2,4,3.3,3.48774919556346])
E2_l  = np.array([2,4,6,8])
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
plt.plot(E2_l,Edot,'o',label = 'gamma =1.36')
plt.grid()
plt.ylabel('$\u03BB/{x_d}$',size=18)
plt.xlabel('effective activation energy',size=18)
plt.legend(loc='upper right',fontsize = 16)
plt.xlim(0,20)
plt.ylim(0,70)
plt.xticks(size=16)
plt.yticks(size=16)
plt.show()