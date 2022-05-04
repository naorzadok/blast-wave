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
gamma = 1.333
q_hat = 26.333
Mw = 0.027
P0 = 101325.0
ood = [1.5,1.6,1.7,1.8,1.9,2.0]
od_num = 1
c = [2,4,6,8,10,12,14,16]
j = 1
gamma_l = np.linspace(1.1,1.7,50,endpoint=True)
print('gamma loop')
d0_l, L_l, W_l, theta_l = [],[],[],[]
#gamma loop
print(dt.datetime.now())
for loop in gamma_l:
    x = blast(lu_et.get(f'eff{c[j]}')[0],c[j],T0,loop,q_hat,Mw,P0,ood[od_num])
    print('!!!length units in centimeters!!!')
    print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
    print( 'od is', x[3] , 'r0 is [cm] -', x[4])
    print( x[4] , x[1] , x[0] , x[2] , x[3])
    d0_l.append(x[4])
    L_l.append(x[0])
    W_l.append(x[1])
    theta_l.append(x[2])
print(dt.datetime.now())
print('finish gamma loop','\n')
plt.plot(gamma_l,d0_l,label = 'd0 vs gamma')
plt.plot(gamma_l,L_l,label = 'L vs gamma')
plt.plot(gamma_l,W_l,label = 'W vs gamma')
plt.grid()
plt.title('gamma',size=14)
plt.ylabel('x, cm',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('gamma',size=14)
plt.legend(loc='upper left',fontsize = 14)
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


#Ea loop
Ea_l = np.linspace(3,6,30,endpoint=True)
print('Ea loop')
d0_l, L_l, W_l, theta_l = [],[],[],[]
print(dt.datetime.now())
for loop in Ea_l:
    x = blast(lu_et.get(f'eff{c[j]}')[0],loop,T0,gamma,q_hat,Mw,P0,ood[od_num])
    print('!!!length units in centimeters!!!')
    print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
    print( 'od is', x[3] , 'r0 is [cm] -', x[4])
    print( x[4] , x[1] , x[0] , x[2] , x[3])
    d0_l.append(x[4])
    L_l.append(x[0])
    W_l.append(x[1])
    theta_l.append(x[2])
print(dt.datetime.now())
print('finish Ea loop','\n')
plt.plot(Ea_l,d0_l,label = 'd0 vs Ea')
plt.plot(Ea_l,L_l,label = 'L vs Ea')
plt.plot(Ea_l,W_l,label = 'W vs Ea')
plt.grid()
plt.title('Ea',size=14)
plt.ylabel('x, cm',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('Ea',size=14)
plt.legend(loc='upper left',fontsize = 14)
plt.savefig('Ea_change.png')
plt.clf()
#plt.show()

plt.plot(Ea_l,theta_l,label = 'theta')
plt.grid()
plt.ylim([0,50])
plt.title('Ea',size=14)
plt.ylabel('degree',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('Ea',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('Ea_theta.png')
plt.clf()
#plt.show()


#od loop
od_l = np.linspace(1.5,2,30,endpoint=True)
print('od loop')
d0_l, L_l, W_l, theta_l = [],[],[],[]
print(dt.datetime.now())
for loop in od_l:
    x = blast(lu_et.get(f'eff{c[j]}')[0],c[j],T0,gamma,q_hat,Mw,P0,loop)
    print('!!!length units in centimeters!!!')
    print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
    print( 'od is', x[3] , 'r0 is [cm] -', x[4])
    print( x[4] , x[1] , x[0] , x[2] , x[3])
    d0_l.append(x[4])
    L_l.append(x[0])
    W_l.append(x[1])
    theta_l.append(x[2])
print(dt.datetime.now())
print('finish od loop','\n')
plt.plot(od_l,d0_l,label = 'd0 vs od')
plt.plot(od_l,L_l,label = 'L vs od')
plt.plot(od_l,W_l,label = 'W vs od')
plt.grid()
plt.title('OverDrive',size=14)
plt.ylabel('x, cm',size=14)
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
plt.ylabel('degree',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('od',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('od_theta.png')
plt.clf()
#plt.show()

#mol loop
mol_l = np.linspace(0.015,0.05,30,endpoint=True)
print('Mw loop')
d0_l, L_l, W_l, theta_l = [],[],[],[]
print(dt.datetime.now())
for loop in mol_l:
    x = blast(lu_et.get(f'eff{c[j]}')[0],c[j],T0,gamma,q_hat,loop,P0,ood[od_num])
    print('!!!length units in centimeters!!!')
    print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
    print( 'od is', x[3] , 'r0 is [cm] -', x[4])
    print( x[4] , x[1] , x[0] , x[2] , x[3])
    d0_l.append(x[4])
    L_l.append(x[0])
    W_l.append(x[1])
    theta_l.append(x[2])
print(dt.datetime.now())
print('finish mol loop','\n')
plt.plot(mol_l,d0_l,label = 'd0 vs mol')
plt.plot(mol_l,L_l,label = 'L vs mol')
plt.plot(mol_l,W_l,label = 'W vs mol')
plt.grid()
plt.title('molecular weight',size=14)
plt.ylabel('x, cm',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('mol, kg/mol',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('mol_change.png')
plt.clf()
#plt.show()

plt.plot(mol_l,theta_l,label = 'theta')
plt.grid()
plt.ylim([0,50])
plt.title('od',size=14)
plt.ylabel('degree',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('mol, kg/mol',size=14)
plt.legend(loc='upper right',fontsize = 14)
plt.savefig('mol_theta.png')
plt.clf()
#plt.show()


#pressure loop
p_l = np.linspace(10000,199000,30,endpoint=True)
print('pressure loop')
d0_l, L_l, W_l, theta_l = [],[],[],[]
print(dt.datetime.now())
for loop in p_l:
    x = blast(lu_et.get(f'eff{c[j]}')[0],c[j],T0,gamma,q_hat,Mw,loop,ood[od_num])
    print('!!!length units in centimeters!!!')
    print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
    print( 'od is', x[3] , 'r0 is [cm] -', x[4])
    print( x[4] , x[1] , x[0] , x[2] , x[3])
    d0_l.append(x[4])
    L_l.append(x[0])
    W_l.append(x[1])
    theta_l.append(x[2])
print(dt.datetime.now())
print('finish pressure loop','\n')
plt.plot(p_l,d0_l,label = 'd0 vs mol')
plt.plot(p_l,L_l,label = 'L vs mol')
plt.plot(p_l,W_l,label = 'W vs mol')
plt.grid()
plt.title('pressure',size=14)
plt.ylabel('x, cm',size=14)
plt.xticks(size=14)
plt.yticks(size=14)
plt.xlabel('pressure, Pa',size=14)
plt.legend(loc='upper left',fontsize = 14)
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

'''
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
