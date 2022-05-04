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

lu_et = { 'eff2' : [190032.6531055641,12.566455531717144], 
        'eff4' :[931273.94,25.13291106343429e0], 
        'eff6' : [4721659.942973013,37.69936659515144 ], 
        'eff7' :[10784134.425409157,43.98259436101001 ],
        'eff8' : [24841031.8513083,50.26582212686858], 
        'eff10' :[135431427.96463335,62.83227765858573], 
        'eff12' : [766173419.0282595,75.39873319030288], 
        'eff14' :[4479235718.399068,87.96518872202002],
        'eff16' :[26993727257.043926,100.53164425373716  ]    }
#    fo.write('M1, p, rho , u1, T1, T_CJ= {:.5E} {:.5E} {:.5E} {:.5E} {:.5E} {:.5E} \n'.format(M1,p,rho,u1,T1,T_cj))

#def blast(AA,eff,T00,gam,qq,moll,P0):
#x = blast(931273.94,4,298.32916984,1.333,26.333,0.027,101325.0,1.75)
#return(length,width,angle,overdrive,kernel- r0)[cm]
#print(x)

T0 = 298.32916984
gamma = 1.333
q_hat = 26.333
Mw = 0.027
P0 = 101325.0
od = 1.6
k = []
print(dt.datetime.now())
c = [2,4,6,7,8,10,12,14,16]
j=0
for i in lu_et:
    x = blast(lu_et.get(i)[0],c[j],T0,gamma,q_hat,Mw,P0,od)
    j+=1
    #print('!!!length units in centimeters!!!')
    #print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
    #print( 'od is', x[3] , 'r0 is [cm] -', x[4])
    #print( x[4] , x[1] , x[0] , x[2] , x[3])
    print('eff ={} ,cell width, cell length'.format(i), x[1], x[0])
    k.append(x[1])
    #print(dt.datetime.now(),i, '\n')
#print(k)
print(dt.datetime.now(), '\n')
'''
results = []
print(datetime.datetime.now())
c = [2,4,6,7,8,10,12,14,16]
j=0
for i in lu_et:
    print(lu_et.get(i)[0],c[j])
    x = blast(lu_et.get(i)[0],c[j],298.32916984,1.333,26.333,0.027,101325.0)
    j+=1
    print('!!!length units in centimeters!!!')
    print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
    print( 'od is', x[3] , 'r0 is [cm] -', x[4], 'ignition tempreture', x[5])
    print( datetime.datetime.now(), i ,'\n')
    results.append(x)
for i in results:
    print(i[4],i[1],i[0],i[2],i[3],i[5])
'''
