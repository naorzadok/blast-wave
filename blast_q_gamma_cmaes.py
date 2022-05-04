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
import optuna


T0 = 300.0
gamma = 1.4
q_hat = 30.0
Mw = 0.036
P0 = 101325.0
ood = 1.6
gamma_l = np.linspace(1.1,1.7,13)
q_l = np.linspace(10,60,17)
c = np.linspace(2,12,11)
aa= 11e4
print(q_l)
print(gamma_l)
print(c)
kkk=np.empty(((len(q_l)*len(gamma_l)),len(c)))#
# kkk= np.empty(shape=(size=(len(gamma_l),len(c))))#.reshape(size=(len(gamma_l),len(c)))#np.empty(len(gamma_l))
#print(kkk)
j=0
for id_q,ps in enumerate(q_l):
    for id_gam,ps_gam in enumerate(gamma_l):
        #print('gamma = ',ps)
        for idx,i in enumerate(c):
            AA=aa*10**(idx*0.75/2)
            x = blast(AA,i,T0,ps_gam,ps,Mw,P0,ood)
            #print('!!!length units in centimeters!!!')
            #print('cell length is [cm] -', x[0], 'cell width is [cm]', x[1],    'cell angle is [deg]', x[2] )
            #print( 'od is', x[3] , 'r0 is [cm] -', x[4])
            #print( x[4] , x[1] , x[0] , x[2] , x[3])
            data1step,Xd_in,eff=zndznd(ps_gam,ps,AA,x[5],1)
            Xd=data1step[0][Xd_in]
            W=x[1]/100
            #print(W,Xd,x[5])
            kkk[j][idx]=W/Xd
        j+=1
print(kkk)

plt.rcParams['font.size'] = '32'
'''
for idd,i in enumerate(kkk):
    x = np.array(c)
    y = np.array(i)
    if idd==round((len(q_l)*len(gamma_l)/2)):
        plt.plot(x,y,label='base case')
    else:
        plt.plot(x,y)
plt.legend() 
plt.grid()
plt.ylabel('W/Xd')
plt.xlabel('effective activation energy')
plt.show()
'''
def gamma_y(a):
    global kk
    cc1=0
    cc2=0
    sumy=0
    kk=np.empty(kkk.shape)
    for jdx,j in enumerate(kkk):
        kk[jdx]=(np.log(kkk[jdx])/np.log(a[2]))*q_l[cc1]**a[0]*gamma_l[cc2]**a[1]
        if (cc2+1)%13==0:
            cc1+=1
            cc2=-1
        cc2+=1
        ##np.log(array) / np.log(base)
    for idx,i in enumerate(kk):
        sumy+=np.sum(abs(kk[0]-i))
    return sumy

x=[2,2,5]

#######cmaes opt
def objective(trial: optuna.Trial):
    x0 = trial.suggest_float("x0", -3.0, 10.0)
    x1 = trial.suggest_float("x1", 0.0, 10.0) #log 
    x2 = trial.suggest_float("x2", 0.01, 10.0)
    x=[x0,x1,x2]
    return gamma_y(x)

if __name__ == "__main__":
    sampler = optuna.samplers.CmaEsSampler(consider_pruned_trials=True)
    study = optuna.create_study(sampler=sampler)
    study.optimize(objective, n_trials=200)
    fig = optuna.visualization.plot_optimization_history(study)
    fig.show()
    print("Best value: {} (params: {})\n".format(study.best_value, study.best_params))
print(study.best_params)
dicty=study.best_params
x00=dicty.get('x0')
x11=dicty.get('x1')
x22=dicty.get('x2')
xx=[x00,x11,x22]
################\

#gam_power =opt.least_squares(gamma_y,x,bounds=((-5,-5,0.01),(10,10,15)))
#gam_power =opt.least_squares(gamma_y,x,method='lm')
#gam_power=opt.fmin_cg(gamma_y,x)
#print(gam_power)
#leftover=gamma_y(gam_power) ###fmin
#gam_power,info, ier, mesg =  opt.fsolve(gamma_y,x,xtol=1.0e-15,full_output=True)
#print(gam_power,info,ier,mesg)
leftover=gamma_y(xx) ###least_squares
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
plt.ylabel('log(W/Xd)*q^a')
plt.xlabel('effective activation energy')
plt.show()
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
