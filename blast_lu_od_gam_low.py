from re import X
from typing import List
import numpy as np
from numpy.core.defchararray import array
from numpy.core.fromnumeric import argmin
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
import datetime as dt
from blast_model import blast

gam_low = { 'eff2' :[165016.21506052298,7.6311988623132905],
          'eff4' :[623163.4686931444,7.6311988623132905*2],
          'eff6' :[2601090.200200781,7.6311988623132905*3],
          'eff8' :[11875856.53597567,7.6311988623132905*4],
          'eff10':[56914457.6900561,7.6311988623132905*5],
          'eff12':[308370047.0914552,7.6311988623132905*6] }
##eff14 :[{'x0': 1700833228.0488868},7.933054597634038e-05], Ea = 53.41839203619303
#'eff4' :[623163.4686931444,7.6311988623132905*2],


T0 = 300.0
gamma = 1.17
q_hat = 50.0
Mw = 2.44412903e-02
P0 = 101325.0
#ood = [1.5,1.6,1.7,1.8,1.9,2.0]
ood=[1.6]
for ps in ood:
    k = []
    #print(dt.datetime.now())
    c = [2,4,6,8,10,12]
    j=0
    for i in gam_low:
        #def blast(AA,eff,T00,gam,qq,moll,P00,odd):
        x = blast(gam_low.get(i)[0],c[j],T0,gamma,q_hat,Mw,P0,ps)
        #print('!!!length units in centimeters!!!')
        #print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
        #print( 'od is', x[3] , 'r0 is [cm] -', x[4])
        #print( x[4] , x[1] , x[0] , x[2] , x[3])
        print('{} ,cell width={}, cell length={},r0={}'.format((i), x[1], x[0],x[4]))
        k.append(x[1])
        #print(dt.datetime.now(),i, '\n')
    print(k)
    #print(dt.datetime.now(), '\n')
'''
        print('!!!length units in centimeters!!!')
        print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
        print( 'od is', x[3] , 'r0 is [cm] -', x[4])
        print( x[4] , x[1] , x[0] , x[2] , x[3])
        k.append(x[1])
        print(dt.datetime.now(),i,c[j], '\n')
        j+=1
    print(k,'\n')
'''

'''
results = []
print(datetime.datetime.now())
c = [2,4,6,7,8,10,12,14,16]
j=0
for i in gam_low:
    print(gam_low.get(i)[0],c[j])
    x = blast(gam_low.get(i)[0],c[j],298.32916984,1.333,26.333,0.027,101325.0)
    j+=1
    print('!!!length units in centimeters!!!')
    print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
    print( 'od is', x[3] , 'r0 is [cm] -', x[4], 'ignition tempreture', x[5])
    print( datetime.datetime.now(), i ,'\n')
    results.append(x)
for i in results:
    print(i[4],i[1],i[0],i[2],i[3],i[5])
'''
