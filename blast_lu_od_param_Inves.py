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
        'eff8' : [24841031.8513083,50.26582212686858], 
        'eff10' :[135431427.96463335,62.83227765858573], 
        'eff12' : [766173419.0282595,75.39873319030288], 
        'eff14' :[4479235718.399068,87.96518872202002],
        'eff16' :[26993727257.043926,100.53164425373716]}
#    fo.write('M1, p, rho , u1, T1, T_CJ= {:.5E} {:.5E} {:.5E} {:.5E} {:.5E} {:.5E} \n'.format(M1,p,rho,u1,T1,T_cj))

#def blast(AA,eff,T00,gam,qq,moll,P0):
#x = blast(931273.94,4,298.32916984,1.333,26.333,0.027,101325.0,1.75)
#return(length,width,angle,overdrive,kernel- r0)[cm]
#print(x)
T0 = 298.32916984
gamma = 1.36
q_hat = 26.333
Mw = 0.027
P0 = 101325.0
ood = [1.5,1.6,1.7,1.8,1.9,2.0]
gamma_l = np.linspace(1.15,1.5,8,endpoint=True)
print(gamma_l)
d0_l, L_l, W_l, theta_l = [],[],[],[]
#gamma loop
for loop in gamma_l:
    k = []
    print(dt.datetime.now())
    c = [2,4,6,8,10,12,14,16]
    j=0
    for i in lu_et:
        x = blast(lu_et.get(i)[0],c[j],T0,loop,q_hat,Mw,P0,ood[1])
        j+=1
        print('!!!length units in centimeters!!!')
        print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
        print( 'od is', x[3] , 'r0 is [cm] -', x[4])
        print( x[4] , x[1] , x[0] , x[2] , x[3])
        k.append(x[1])
        d0_l.append(x[4])
        L_l.append(x[0])
        W_l.append(x[1])
        theta_l.append(x[2])
        print(dt.datetime.now(),i, '\n')
        #print(k)
    plt.plot(gamma_l,d0_l,label = f'd0 vs gamma for gamma{loop}')
    plt.plot(gamma_l,L_l,label = f'L vs gamma for gamma{loop}')
    plt.plot(gamma_l,W_l,label = f'W vs gamma for gamma{loop}')
print('finish gamma loop')
d0_l, L_l, W_l, theta_l = [],[],[],[]
#plt.plot(gamma_l,d0_l,label = 'd0 vs gamma')
#plt.plot(gamma_l,L_l,label = 'L vs gamma')
#plt.plot(gamma_l,W_l,label = 'W vs gamma')
plt.grid()
plt.title('gamma',size=14)
plt.ylabel('[m]',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('gamma',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('gamma_change.png')
plt.clf()
#plt.show()

plt.plot(gamma_l,theta_l,label = 'theta')
plt.grid()
plt.ylim([0,50])
plt.title('gamma',size=14)
plt.ylabel('[degree]',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('gamma',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('gamma_theta.png')
plt.clf()
#plt.show()

d0_l, L_l, W_l, theta_l = [],[],[],[]
#Ea loop
Ea_l = np.linspace(20,40,8,endpoint=True)
print(Ea_l)
for i in Ea_l:
    d0,L,W,theta = main_mod(1.333,i,1.794,0.027,1.0)
    d0_l.append(d0)
    L_l.append(L)
    W_l.append(W)
    theta_l.append(theta)
print('finish Ea loop')

plt.plot(Ea_l,d0_l,label = 'd0 vs Ea')
plt.plot(Ea_l,L_l,label = 'L vs Ea')
plt.plot(Ea_l,W_l,label = 'W vs Ea')
plt.grid()
plt.title('Ea',size=14)
plt.ylabel('[m]',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('Ea',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('Ea_change.png')
plt.clf()
#plt.show()

plt.plot(Ea_l,theta_l,label = 'theta')
plt.grid()
plt.ylim([0,50])
plt.title('Ea',size=14)
plt.ylabel('[degree]',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('Ea',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('Ea_theta.png')
plt.clf()
#plt.show()


d0_l, L_l, W_l, theta_l = [],[],[],[]
#od loop
od_l = np.linspace(1.5,2,8,endpoint=True)
for i in od_l:
    d0,L,W,theta = main_mod(1.333,25.13291106343429,i,0.027,1.0)
    d0_l.append(d0)
    L_l.append(L)
    W_l.append(W)
    theta_l.append(theta)
print('finish od loop')

plt.plot(od_l,d0_l,label = 'd0 vs od')
plt.plot(od_l,L_l,label = 'L vs od')
plt.plot(od_l,W_l,label = 'W vs od')
plt.grid()
plt.title('OverDrive',size=14)
plt.ylabel('[m]',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('od',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('od_change.png')
plt.clf()
#plt.show()

plt.plot(od_l,theta_l,label = 'theta')
plt.grid()
plt.ylim([0,50])
plt.title('od',size=14)
plt.ylabel('[degree]',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('od',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('od_theta.png')
plt.clf()
#plt.show()

d0_l, L_l, W_l, theta_l = [],[],[],[]
#mol loop
mol_l = np.linspace(0.017,0.04,8,endpoint=True)
for i in mol_l:
    d0,L,W,theta = main_mod(1.333,25.13291106343429,1.794,i,1.0)
    d0_l.append(d0)
    L_l.append(L)
    W_l.append(W)
    theta_l.append(theta)
print('finish mol loop')

plt.plot(mol_l,d0_l,label = 'd0 vs mol')
plt.plot(mol_l,L_l,label = 'L vs mol')
plt.plot(mol_l,W_l,label = 'W vs mol')
plt.grid()
plt.title('molecular weight',size=14)
plt.ylabel('[m]',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('mol',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('mol_change.png')
plt.clf()
#plt.show()

plt.plot(mol_l,theta_l,label = 'theta')
plt.grid()
plt.ylim([0,50])
plt.title('od',size=14)
plt.ylabel('[degree]',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('mol',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('mol_theta.png')
plt.clf()
#plt.show()

d0_l, L_l, W_l, theta_l = [],[],[],[]
#pressure loop
p_l = np.linspace(0.5,1.5,11,endpoint=True)
for i in p_l:
    d0,L,W,theta = main_mod(1.333,25.13291106343429,1.794,0.027,i)
    d0_l.append(d0)
    L_l.append(L)
    W_l.append(W)
    theta_l.append(theta)
print('finish pre loop')
plt.plot(p_l,d0_l,label = 'd0 vs mol')
plt.plot(p_l,L_l,label = 'L vs mol')
plt.plot(p_l,W_l,label = 'W vs mol')
plt.grid()
plt.title('molecular weight',size=14)
plt.ylabel('[m]',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('mol',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('p_change.png')
plt.clf()
#plt.show()

plt.plot(p_l,theta_l,label = 'theta')
plt.grid()
plt.ylim([0,50])
plt.title('od',size=14)
plt.ylabel('[degree]',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('mol',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('p_theta.png')
plt.clf()
#plt.show()

T0 = 298.32916984
gamma = 1.36
q_hat = 26.333
Mw = 0.027
P0 = 101325.0
ood = [1.5,1.6,1.7,1.8,1.9,2.0]
for ps in ood:
    k = []
    print(dt.datetime.now())
    c = [2,4,6,7,8,10,12,14,16]
    j=0
    for i in lu_et:
        x = blast(lu_et.get(i)[0],c[j],T0,gamma,q_hat,Mw,P0,ps)
        j+=1
        print('!!!length units in centimeters!!!')
        print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
        print( 'od is', x[3] , 'r0 is [cm] -', x[4])
        print( x[4] , x[1] , x[0] , x[2] , x[3])
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
