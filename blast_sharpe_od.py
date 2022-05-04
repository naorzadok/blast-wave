from os import write
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
from blast_model_sharpe import blast

sharpe_et = { 'eff2' :[480.8413988062197,13.333333333333332],
              'eff3' :[1154.8916494693397,20.0],
              'eff4' :[3175.2929102699927, 26.666666666666664],
              'eff6' :[17869.308484519064, 40.0],
              'eff8' :[100127.9651722467, 53.33333333333333],
              'eff10' :[581413.7985920393, 66.66666666666666]}
#    fo.write('M1, p, rho , u1, T1, T_CJ= {:.5E} {:.5E} {:.5E} {:.5E} {:.5E} {:.5E} \n'.format(M1,p,rho,u1,T1,T_cj))

#def blast(AA,eff,T00,gam,qq,moll,P0):
#x = blast(931273.94,4,298.32916984,1.333,26.333,0.027,101325.0,1.75)
#return(length,width,angle,overdrive,kernel- r0)[cm]
#print(x)
f = open('sharpe_odres.csv','w')
T0 = 300.0
gamma = 1.54
q_hat = 11.5
Mw = 0.031568
P0 = 101325.0
ood = [1.6]
for ps in ood:
    k = []
    print(dt.datetime.now())
    c = [2,3,4,6,8,10]
    j=0
    f.write(f' od={ps}, \n')
    for i in sharpe_et:
        x = blast(sharpe_et.get(i)[0],c[j],T0,gamma,q_hat,Mw,P0,ps)
        j+=1
        print('!!!length units in centimeters!!!')
        print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
        print( 'od is', x[3] , 'r0 is [cm] -', x[4] ,'\n')
        #print( x[4] , x[1] , x[0] , x[2] , x[3])
        f.write(f'{x[4]} , {x[1]} , {x[0]} , {x[2]} , {x[3]} ,\n')
        k.append(x[1])
        print(dt.datetime.now(),i, '\n')
    print(k)
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
