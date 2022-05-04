from cProfile import label
from re import X
#from turtle import color
#from tkinter import W
from typing import List
from matplotlib.transforms import BboxBase
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
from fractions import Fraction
from IPython.display import display, Math, Latex

# units in centimeters
########sharpe###########
#Ea3    = np.array([9.631248145576679])
Ea3 = np.array([7.286722469170652, 10.303713848332231, 14.272126621536309, 18.934320956368634, 23.893216946469305])
E_l_sh  = np.array([3])
E_res = np.array([80/7])
E_l_sharpe  = np.array([2,4,6,8,10])
Edot    = E_res/1.0
##########################

####high gamma###
#W16_l_high = np.array([1.3142224876525281, 2.2531968317912887, 3.1461967177948926, 4.112009273030238])
W16_l_high = np.array([1.3145005934249303, 2.2533265176600508, 3.1473177700386983, 4.110565750507097, 5.1971036254122955, 6.420863583536697])
E_res = np.array([2.1333333,3.76470588235294,3.3])/0.175
E_l  = np.array([2,4,6,8,10,12])
E2_l  = np.array([2,4,6])
x16_ax_high  = W16_l_high/0.175
#################

####high q###
W16_l_high_q=np.array([2.0152145274116444, 3.87175500990289, 6.017460302417044, 7.2457688704033165, 8.600489665377692, 11.694812098770406, 15.256711033741846, 19.232617971300748, 23.509030630275348])#
#[2.2283549929063122, 4.436693421408165, 7.072866967085392, 10.276651460191752, 14.12380248850615, 18.569446147381107])
E_l_high_q = np.array([2,4,6,8,10,12])
x16_ax_high_q  = W16_l_high_q/0.175
E_high_q_res=np.array([2.56,32/8,3.2])/0.175
E2_l_high_q  = np.array([2,4,6])
#####################

#########lu et al.###########
W16_l = np.array([2.0099448159652655, 3.8616304616834367, 6.0017247326066645, 8.577999547773649, 11.664230325132511, 15.216815094569808, 19.182325129368397, 23.447555070072138])
x16_ax = W16_l/0.175
E_l_lu  = np.array([2,4,6,8,10,12,14,16])
Edot_lu = np.array([17.274618317056497,23.728294177732373,34.111531502363086,28.013687149637473,27.546570548130653,28.01766980360358,22.917351334547227,17.523004124482284])
############################

##########q_up#############
W16_l_q_up = np.array([2.0079095789135577, 3.8971280594764726, 6.1017728461465195, 8.771413803600637, 11.975003922447678, 15.682001474450658, 19.811107257088878, 24.248338913680286])
E_l_q_up  = np.array([2,4,6,8,10,12,14,16])
E_res_q_up = np.array([4.0])
E2_l_q_up  = np.array([4.0])
x_ind = 0.175
#x15_ax  = W15_l/x_ind
x16_ax_q_up  = W16_l/x_ind
Edot_q_up    = E_res_q_up/x_ind

'''
filename='error_barz.csv'
fo = open(filename,'w')

fo.write(f'Ea-energy, \u03B3, Q=(qMw/RT), \u03BB-simulation cm, \u03BB-Our model cm, error % \n')
fo.write(f'4, 1.333, 26.333, {round(23.728294177732373*0.175,2)}, 3.8616304616834367, {round(np.abs((23.728294177732373*0.175-3.8616304616834367)/(23.728294177732373*0.175))*100,2)} \n')
fo.write(f'6, 1.333, 26.333, {round(34.111531502363086*0.175,2)}, 6.0017247326066645, {round(np.abs((34.111531502363086*0.175-6.0017247326066645)/(34.111531502363086*0.175))*100,2)}\n')
fo.write(f'6, 1.5, 26.333, {round(3.1473177700386983,2)}, 3.30, {round(np.abs((3.1473177700386983-3.3)/(3.1473177700386983))*100,2)}\n')
fo.write(f'4, 1.333, 30.0, {round(3.8971280594764726,2)}, 4.0, {round(np.abs((3.8971280594764726-4.0)/(3.8971280594764726))*100,2)}\n')
fo.write(f'3, 1.54, 11.5, {round(9.631248145576679,2)}, {round(80/7,2)}, {round(np.abs((9.631248145576679-80/7)/(80/7))*100,2)}\n')
fo.close()
'''
#for k,val in enumerate(X):
#    for j in X[k]:
#        fo.write(f'{j} ,')
#    fo.write(f' , {F[k,0]} , {F[k,1]}  \n')



plt.rcParams['legend.title_fontsize'] = 20#'large'
plt.rcParams['lines.markersize'] = 14
plt.rcParams['lines.linewidth'] = 5
plt.rcParams['legend.fontsize'] = 20
plt.rcParams['lines.markeredgewidth'] = 3.0
#plt.plot(E_l,x15_ax,label = 'overdrive 1.5')

################lu
j1, =plt.plot(E_l_lu,x16_ax,'b--')#, label='\u03B3=1.333,Q=26.333')
j2, =plt.plot(E_l_lu,Edot_lu,'bo',markeredgecolor='black')#,label = 'Lu et al')
plt.errorbar(E_l_lu,Edot_lu,yerr=2.5,fmt = 'o',ecolor='b',color='b',markeredgecolor='black',elinewidth = 4, capsize=5)
plt.gca().add_artist(plt.legend([j1,j2],['Our model','simulation'],title='Lu et al. :$\gamma$=1.333; Q=26.333',loc='upper left'))#,bbox_to_anchor=(1.0,0.75)))#,fontsize=18))
################

############high_gamma
l1, =plt.plot(E_l,x16_ax_high,'g--',markersize=20)#,label = 'Our model model')
l2, =plt.plot(E2_l,E_res,'go',markeredgecolor='black')#3,label = 'simulation')
plt.errorbar(E2_l,E_res,yerr=2.5,fmt = 'o',ecolor='g',color='g',markeredgecolor='black',elinewidth = 4, capsize=5)
plt.gca().add_artist(plt.legend([l1,l2],['Our model','simulation'],title='high $\gamma$ =1.5 ; Q=26.333',loc=0,bbox_to_anchor=(1.0,1.0)))#,fontsize=18))
############high_q
f1, =plt.plot(E_l_high_q,x16_ax_high_q,'y--',markersize=20)#,label = 'Our model model')
f2, =plt.plot(E2_l_high_q,E_high_q_res,'yo',markeredgecolor='black')#3,label = 'simulation')
plt.errorbar(E2_l_high_q,E_high_q_res,yerr=2.5,fmt = 'o',ecolor='y',color='y',markeredgecolor='black',elinewidth = 4, capsize=5)
plt.gca().add_artist(plt.legend([f1,f2],['Our model','simulation'],title='high q, $\gamma$ =1.333 ; Q=50.0',loc=1,bbox_to_anchor=(1.0,0.8)))
################q_up
t1, =plt.plot(E_l_q_up,x16_ax_q_up,'r-.')#, label='\u03B3=1.333,Q=26.333')
t2, =plt.plot(E2_l_q_up,Edot_q_up,'ro',markeredgecolor='black')#,label = 'Lu et al')
plt.errorbar(E2_l_q_up,Edot_q_up,yerr=2.5,ecolor='r',elinewidth = 4, capsize=5)
plt.gca().add_artist(plt.legend([t1,t2],['Our model','simulation'],title='$Q_{up}$ :$\gamma$=1.333 ; Q=30.0',loc='lower right'))#,bbox_to_anchor=(1.0,0.5)))#,fontsize=18))
######################sharpe
k1, =plt.plot(E_l_sharpe,Ea3,'m--')#,label = 'sharpe model')
k2, =plt.plot(E_l_sh,Edot,'mo',markeredgecolor='black')#,label = 'sharpe simulation')
#plt.errorbar(E_l_sh,Edot,yerr=1/(32*0.175),ecolor='m',elinewidth = 4, capsize=5)
plt.gca().add_artist(plt.legend([k1,k2],['Our model','simulation'],title='Sharpe et al.-$\gamma$=1.54 ; Q=11.5',loc='lower center'))#,bbox_to_anchor=(0.75,0.2)))#,fontsize=18))
#####################


plt.grid()
plt.ylabel('$\u03BB/{x_d}$',size=22)
plt.xlabel('effective activation energy',size=22)
#plt.legend(loc='upper right',fontsize = 16)
plt.xlim(0,16)
plt.ylim(0,50)
plt.xticks(size=22)
plt.yticks(size=22)
plt.show()