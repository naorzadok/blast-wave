from re import X
from typing import List
import numpy as np
from numpy import size
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
from blast_model_eff import blast
from znd_fun_noWrite import zndznd


T0 = 300.0
gamma = 1.4
q_hat = 30.0
Mw = 0.036
P0 = 101325.0
ood = 1.6
gamma_l = np.linspace(1.1,1.7,7)
q_l = np.linspace(10,60,6)
c = np.linspace(2,12,11)
aa= 11e4
print(gamma_l)
print(c)
'''
kkk=np.empty((len(gamma_l),len(c)))#
# kkk= np.empty(shape=(size=(len(gamma_l),len(c))))#.reshape(size=(len(gamma_l),len(c)))#np.empty(len(gamma_l))
#print(kkk)
j=0
for id,ps in enumerate(gamma_l):
    k = []
    #print(dt.datetime.now())
    #print('gamma = ',ps)
    for idx,i in enumerate(c):
        AA=aa*10**(idx*0.75/2)
        x = blast(AA,i,T0,ps,q_hat,Mw,P0,ood)
        #print('!!!length units in centimeters!!!')
        #print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
        #print( 'od is', x[3] , 'r0 is [cm] -', x[4])
        #print( x[4] , x[1] , x[0] , x[2] , x[3])
        data1step,Xd_in,eff=zndznd(ps,q_hat,AA,x[5],1)
        Xd=data1step[0][Xd_in]
        W=x[1]/100
        #print(W,Xd,x[5])
        kkk[id][idx]=W/Xd
print(kkk)
'''
def fitt(a):
    sumy=0
    kk=np.empty((len(gamma_l),31))
    for idd,i in enumerate(kkk):
        x = np.array(np.power(c,a[0]))*gamma_l[idd]**a[1]
        y = np.array(i)
        kk[idd]=np.polyfit(x, y, 30)
    for idx,i in enumerate(kk):
        sumy+=np.sum(abs(kk[0]-i))
    print(kk)
    return sumy
    
plt.rcParams['font.size'] = '32'
kkk=np.load('gamma_l_file.npy')
a=-0.75#-0.5
b=3#-0.5
x=[a,b]
xx=x
#xx=opt.least_squares(fitt,x)#,max_nfev=1)
print(xx)
#xx=xx.x
for idd,i in enumerate(kkk):
    x = np.array(np.power(c,xx[0]))*gamma_l[idd]**xx[1]
    #x = np.array(np.power(c,a))*gamma_l[idd]**b
    y = np.array(i)
    plt.plot(x,y, label='gamma={}'.format(1.0 +0.1*(idd+1)))
    #if idd==0:
    #    z=np.poly1d(np.polyfit(x, y, 6))
    #    plt.plot(x,z(x),label='polyfit')
plt.legend() 
#plt.ylim(0,200)
#plt.xlim(0,4)
plt.grid()
plt.ylabel('W/Xd')
plt.xlabel('\u03B5^{}*gamma^{}'.format(xx[0],xx[1]))
#plt.xlabel('\u03B5^{}*gamma^{}'.format(a,b))
plt.show()

for idd,i in enumerate(kkk):
    x = np.array(c)
    y = np.array(i)
    plt.plot(x,y, label='gamma={}'.format(1.0 +0.1*(idd+1)))
    #if idd==6:
    #    plt.plot(x,y,label='base case')
    #else:
    #    plt.plot(x,y)
plt.legend() 
plt.grid()
plt.ylabel('W/Xd')
plt.xlabel('\u03B5')
plt.show()





'''
print(kkk)
plt.rcParams['font.size'] = '32'
for idd,i in enumerate(kkk):
    x = np.array(c)
    y = np.array(i)
    plt.plot(x,y)
    #if idd==6:
    #    plt.plot(x,y,label='base case')
    #else:
    #    plt.plot(x,y)
plt.legend() 
plt.grid()
plt.ylabel('W/Xd')
plt.xlabel('effective activation energy')
plt.show()

def gamma_y(a):
    global kk
    sumy=0
    kk=np.empty(kkk.shape)
    for jdx,j in enumerate(kkk):
        kk[jdx]=(np.log(kkk[jdx])/np.log(a[1]))*gamma_l[jdx]**a[0]
        ##np.log(array) / np.log(base)
    for idx,i in enumerate(kk):
        sumy+=np.sum(abs(kk[0]-i))
    return [sumy,sumy]

x=[2,4]
#gam_power =opt.root_scalar(gamma_y,method='secant',x0=0.004,x1=0.002)
gam_power,info, ier, mesg =  opt.fsolve(gamma_y,x,xtol=1.0e-15,full_output=True)
print(gam_power,info,ier,mesg)
leftover=gamma_y(gam_power)
print(kk,leftover)
for idd,i in enumerate(kk):
    x = np.array(c)
    y = np.array(i)
    if idd==6:
        plt.plot(x,y,label='base case')
    else:
        plt.plot(x,y)
plt.legend() 
plt.grid()
plt.ylabel('log(W/Xd)*gamma^a')
plt.xlabel('effective activation energy')
plt.show()
'''
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
