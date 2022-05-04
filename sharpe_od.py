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

sharpe_et = { 'eff3' : [1158.6349259643584,20.0]}
#    fo.write('M1, p, rho , u1, T1, T_CJ= {:.5E} {:.5E} {:.5E} {:.5E} {:.5E} {:.5E} \n'.format(M1,p,rho,u1,T1,T_cj))

#def blast(AA,eff,T00,gam,qq,moll,P0):
#x = blast(931273.94,4,298.32916984,1.333,26.333,0.027,101325.0,1.75)
#return(length,width,angle,overdrive,kernel- r0)[cm]
#print(x)

T0 = 300.0
gamma = 1.54
q_hat = 11.5
Mw = 0.031568
P0 = 101325.0
od = 1.6
L_l,L_l3 = [] , []
print(dt.datetime.now())
c = [3]
j=0
od_l = np.linspace(1.5,2.0,6)
k = []
for od in od_l:
    x = blast(sharpe_et.get('eff3')[0],3,T0,gamma,q_hat,Mw,P0,od)
    j+=1
    print('!!!length units in centimeters!!!')
    print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
    print( 'od is', x[3] , 'r0 is [cm] -', x[4])
    print( x[4] , x[1] , x[0] , x[2] , x[3],'\n')
    L_l.append(x[1])
    L_l3.append(x[0])
    k.append(x[1])
print(k,'\n')

L_l2 = np.array(L_l)
ans_l = np.full_like(od_l,80/7)
ans2_l = np.full_like(od_l,11.4)
plt.figure(figsize=(7, 5))
plt.plot(od_l,L_l2,label = 'model results')
plt.plot(od_l,ans2_l,label= 'simulation results')
plt.grid()
plt.xticks(size=14)
plt.yticks(size=14)
#plt.ylim(14,16)
#plt.xlim(1.65,1.7)
plt.ylabel('Cell width, Xd',size=14)
plt.xlabel('overdrive factor',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.show()
#plt.savefig('sharpe_od_zoom.png')


L_l2 = np.array(L_l3)
ans_l = np.full_like(od_l,80/7)
ans2_l = np.full_like(od_l,15)
plt.figure(figsize=(7, 5))
plt.plot(od_l,L_l2,label = 'model results')
plt.plot(od_l,ans2_l,label= 'simulation results')
plt.grid()
plt.xticks(size=14)
plt.yticks(size=14)
#plt.ylim(14,16)
#plt.xlim(1.65,1.7)
plt.ylabel('Cell length, Xd',size=14)
plt.xlabel('overdrive factor',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.show()

plt.figure(figsize=(7, 5))
plt.plot(od_l,L_l,label = 'model results')
plt.plot(od_l,ans2_l,label= 'simulation results')
plt.grid()
plt.xticks(size=14)
plt.yticks(size=14)
plt.ylim(5.0,10.0)
plt.xlim(1.5,2.0)
plt.ylabel('r, cm',size=14)
plt.xlabel('overdrive factor',size=14)
plt.legend(loc='upper right',fontsize = 14)
#plt.show()
#plt.savefig('for_conf_cm.png')
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
