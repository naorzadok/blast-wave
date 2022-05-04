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
#W15_l  = np.array([2.5835715319329986, 5.318868383942599, 8.73771337636946, 13.02078146921781, 18.200332912205386, 24.24368818491984, 30.86856619161073, 37.939908003139216])
W16_l = np.array([2.0079095789135577, 3.8971280594764726, 6.1017728461465195, 8.771413803600637, 11.975003922447678, 15.682001474450658, 19.811107257088878, 24.248338913680286])
#W17_l  = np.array([2.217742976906275, 4.347116258044583, 6.948328511854268, 10.204871959967027, 14.152036620084106, 18.760968253586675, 23.80562788688553, 29.175104984922708])
#W18_l  = np.array([2.082580875487875, 4.009014206889886, 6.345757775137801, 9.270167841158038, 12.812472896594448, 16.941369123959905, 21.447738880026886, 26.2282351526996])
#W19_l  = np.array([1.9680317716212439, 3.729894280411787, 5.853278082900882, 8.506683930318415, 11.71357707545841, 15.439148199049715, 19.48749764061225, 23.761308201392534])
#W20_l  = np.array([1.8682280039516572, 3.4914463267024622, 5.434280314796634, 7.854551052553695, 10.767986757434615, 14.1356480784723, 17.77252899585113, 21.586431332658176])
E_l  = np.array([2,4,6,8,10,12,14,16])
E_res = np.array([4])
E2_l  = np.array([4])
x_ind = 0.175
#x15_ax  = W15_l/x_ind
x16_ax  = W16_l/x_ind
#x17_ax  = W17_l/x_ind
#x18_ax  = W18_l/x_ind
#x19_ax  = W19_l/x_ind
#x20_ax  = W20_l/x_ind
Edot    = E_res/x_ind

eff  = np.linspace(2,20,19)
res  = f_el(eff)


#plt.plot(E_l,x15_ax,label = 'overdrive 1.5')
plt.plot(E_l,x16_ax,label = 'overdrive 1.6')
#plt.plot(E_l,x17_ax,label = 'overdrive 1.7')
#plt.plot(E_l,x18_ax,label = 'overdrive 1.8')
#plt.plot(E_l,x19_ax,label = 'overdrive 1.9')
#plt.plot(E_l,x20_ax,label = 'overdrive 2.0')
#plt.plot(eff,res,label = 'Lu et al')
plt.plot(E2_l,Edot,'o',label = 'q=30.0')
plt.grid()
plt.ylabel('$\u03BB/{x_d}$',size=18)
plt.xlabel('effective activation energy',size=18)
plt.legend(loc='upper right',fontsize = 16)
plt.xlim(0,20)
plt.ylim(0,70)
plt.xticks(size=16)
plt.yticks(size=16)
plt.show()