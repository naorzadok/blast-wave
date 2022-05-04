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

def r(t): #call only after defining t_star
    if t>t_star:
        return np.sqrt(2*D0*np.sqrt(r_star**2+r0**2*(od**2-1))*(t-t_star)+r_star**2)#np.sqrt(2*np.sqrt(2)*r0*D0*od*(t-t_star)+r_star**2)
    else:
        #print('2')
        return np.sqrt((D0*t+r0*od)**2-(r0**2*(od**2-1)))#0.5*np.sqrt((2*D0*t+d0*od)**2-(d0**2*(od**2-1)))

def D(t):
    if t>t_star:
        return (D0/r(t))*np.sqrt(r_star**2+(r0**2*(od**2-1)))
    else:
        #print('1')
        return (D0/r(t))*np.sqrt(r(t)**2+r0**2*(od**2-1))  #np.sqrt((D*d0)**2+4*D0**2*(r(t)**2-r0**2))/(2*r(t))

def D_der(t):
    A = np.sqrt(2*np.sqrt(2)*r0*D0*od)
    return -1/2*(0.5*A*(A(t-t_star)+r_star**2)**(-3/2)*A)

def D_integ(t1,tau):
    return (1/tau)*(quad(D,t1,t1+tau))[0]

def M(t):
    return D(t)/C0

def rho_vn(t):
    return (((gamma+1)*M(t)**2)/(2+(gamma-1)*M(t)**2)*rho0)
def P_vn(t):
    return ((1+(2*gamma)/(gamma+1)*(M(t)**2-1)))*P0
def T_vn(t):
    return (P_vn(t)/(rho_vn(t)*RR))

def tau(t):
    return (1/(A*rho_vn(t))*(np.exp(Ea/(Rr*T_vn(t))))) -  (1/(A*rho_vn0)*(np.exp(Ea/(Rr*T_vn0)))) #reduction of induction time from the temp front
def u1(t):
    return D(t)*((2+(gamma-1)*M(t)**2)/((gamma+1)*M(t)**2))  ###after shock velocity(relative to shock)
def x_ind(t):
    #return tau(t)*D(t)*((gamma-1)/(gamma+1))
    return tau(t)*(D_integ(t,tau(t)))*(1-(gamma-1)/(gamma+1))# integrating on D an has an effect!!!!!!!!1
#def x_ind2(t):
#    return tau(t)*(D_integ(t,tau(t)))*((2+(gamma-1)*M(t)**2)/((gamma+1)*M(t)**2)) #tau(t)*D(t)*(1-(2+(gamma-1)*M(t)**2)/((gamma+1)*M(t)**2)) ##works(after shock absolute velocity) from density ratio

def x_ind2(t):
    return tau(t)*(D_integ(t,tau(t)))*(1-(2+(gamma-1)*M(t)**2)/((gamma+1)*M(t)**2))# uses absolute velocity

def x_ind3(t):
    DD = lambda t1: (D(t1))*(1-(2+(gamma-1)*M(t1)**2)/((gamma+1)*M(t1)**2))
    return (quad(DD,t,t+tau(t)))[0]# uses absolute velocity

def x_ind_rel(t):
    DD = lambda t1: (D(t1))*((2+(gamma-1)*M(t1)**2)/((gamma+1)*M(t1)**2))
    return (quad(DD,t,t+tau(t)))[0]# uses relative velocity
    

def t_final(t): #calculate r at t[0,t*]
    return D0*t+r0-np.sqrt(2*D0*np.sqrt(r_star**2+r0**2*(od**2-1))*(t-t_star)+r_star**2)

def find_t_star(t_star):
    r_star    = np.sqrt((D0*t_star+r0*od)**2-(r0**2*(od**2-1)))
    D_star    = (D0/r_star)*np.sqrt(r_star**2+r0**2*(od**2-1))
    M_star    = D_star/C0
    rho_star  = (((gamma+1)*M_star**2)/(2+(gamma-1)*M_star**2)*rho0)
    P_star    = ((1+(2*gamma)/(gamma+1)*(M_star**2-1)))*P0
    T_star    = (P_star/(rho_star*RR))
    tau       = (1/(A*rho_star)*(np.exp(Ea/(Rr*T_star)))) -  (1/(A*rho_vn0)*(np.exp(Ea/(Rr*T_vn0))))
    t_f       = lambda t: D0*t+r0-np.sqrt(2*D0*np.sqrt(r_star**2+r0**2*(od**2-1))*(t-t_star)+r_star**2)
    #print(np.abs((t_star+tau)-opt.fsolve(t_f,2*t_star)),t_star)
    #print(t_star)
    return ((t_star+tau)-opt.fsolve(t_f,t_ft*t_star))

########################################
def geo_mod(d0): #calculate r at t[0,t*]
    '''
    ############################################
    global r0,t_star, r_star
    r0        = d0/2
    ############################################
    t_star    = opt.root(find_t_star,30*r0/D0).x#t_star    = ((np.sqrt(2)-1))*(od*r0/D0)opt.bisect(find_t_star,0.1*r0/D0,10*r0/D0)#opt.minimize(find_t_star,30*r0/D0,method='nelder-mead').x#
    r_star    = r(t_star)#(r0)*np.sqrt((1+od**2))
    t_end      = opt.fsolve(t_final,1.3*t_star)

    time_jump  = t_star/1000
    t_arr      = np.arange(0, t_end+time_jump,time_jump) # t_arr      = np.arange(0,t_end+time_jump,time_jump)
    #############################################
    r_l       = np.fromiter(map(r, t_arr), dtype=np.float)
    #D_l       = np.fromiter(map(D, t_arr), dtype=np.float)
    #M_l       = np.fromiter(map(M, t_arr),dtype=np.float)
    #rho_l     = np.fromiter(map(rho_vn, t_arr),dtype=np.float)
    #P_l       = np.fromiter(map(P_vn, t_arr), dtype=np.float)
    #T_l       = np.fromiter(map(T_vn, t_arr), dtype=np.float)
    tau_l     = np.fromiter(map(tau, t_arr), dtype=np.float)
    #u1_l      = np.fromiter(map(u1, t_arr), dtype=np.float)
    #x_ind_l   = np.fromiter(map(x_ind, t_arr), dtype=np.float)
    #x_ind2_l  = np.fromiter(map(x_ind2, t_arr), np.float)
    x_ind3_l  = np.fromiter(map(x_ind3, t_arr), np.float)

    t_temp     = tau_l+t_arr
    temp_front = r_l +x_ind3_l
    ix = np.abs(t_temp-t_arr[-1]).argmin()
    temp_front = temp_front[:ix+1]
    t_temp =t_temp[:ix+1]
    print(r0,r0 - (r_l[-1]-temp_front[-1]))
    #print(ix,np.abs((t_arr-t_star)).argmin(), t_star+tau_l[ix+1]-t_arr[-1],t_arr[ix]+tau_l[ix+1]-t_arr[-1])
    return r0 - (r_l[-1]-temp_front[-1])#np.abs(r0-x_ind3[inx_stop])#r0 - (r_arr_f[-1]-temp_front[ix])
    '''
    ############################################
    global r0,t_star, r_star
    r0        = d0/2
    ############################################
    t_star    = opt.root(find_t_star,x0_star*r0/D0).x#t_star    = ((np.sqrt(2)-1))*(od*r0/D0)opt.bisect(find_t_star,0.1*r0/D0,10*r0/D0)#opt.minimize(find_t_star,30*r0/D0,method='nelder-mead').x#
    #t_star    = opt.minimize_scalar(find_t_star,bounds=(1e-8,2),method='brent').x
    r_star    = r(t_star)#(r0)*np.sqrt((1+od**2))
    #print('imma out')
    t_end      = opt.fsolve(t_final,t_ft*t_star)
    if (t_end<1e-12) or (t_star<0):
        return 100
    time_jump  = t_star/1000
    t_arr      = np.arange(0, t_end+time_jump,time_jump) # t_arr      = np.arange(0,t_end+time_jump,time_jump)
    #############################################
    r_l       = np.fromiter(map(r, t_arr), dtype=np.float)
    #D_l       = np.fromiter(map(D, t_arr), dtype=np.float)
    #M_l       = np.fromiter(map(M, t_arr),dtype=np.float)
    #rho_l     = np.fromiter(map(rho_vn, t_arr),dtype=np.float)
    #P_l       = np.fromiter(map(P_vn, t_arr), dtype=np.float)
    #T_l       = np.fromiter(map(T_vn, t_arr), dtype=np.float)
    #tau_l     = np.fromiter(map(tau, t_arr), dtype=np.float)
    #u1_l      = np.fromiter(map(u1, t_arr), dtype=np.float)
    #x_ind_l   = np.fromiter(map(x_ind, t_arr), dtype=np.float)
    #x_ind2_l  = np.fromiter(map(x_ind2, t_arr), np.float)
    #x_ind3_l  = np.fromiter(map(x_ind3, t_arr), np.float)
    x_ind_r = x_ind_rel(t_star)
    '''
    t_temp     = tau_l+t_arr
    temp_front = r_l +x_ind3_l
    ix = np.abs(t_temp-t_arr[-1]).argmin()
    temp_front = temp_front[:ix+1]
    t_temp =t_temp[:ix+1]
    '''
    #print(r0, r0-x_ind_r)#r0 - (r_l[-1]-temp_front[-1]))
    #print(ix,np.abs((t_arr-t_star)).argmin(), t_star+tau_l[ix+1]-t_arr[-1],t_arr[ix]+tau_l[ix+1]-t_arr[-1])
    return  r0-x_ind_r#(r0 - (r_l[-1]-temp_front[-1]))#np.abs(r0-x_ind3[inx_stop])#r0 - (r_arr_f[-1]-temp_front[ix])

def startup(x,gamm,mo,qu,TT):
    T0 = TT
    T1 = x[0]
    T2 = x[1]
    gamma = gamm
    Dcj = x[2]
    RR = 8.3144087/mo
    qq = qu*RR*T0
    Cp = (gamma*RR)/(gamma-1)
    eq1 = T1/T0 -1-((2*(gamma-1))/(gamma+1)**2)*((gamma*((Dcj**2/(gamma*RR*T0))-1)+1-(gamma*RR*T0)/(Dcj**2)))
    eq2 = T2-T0-qq/Cp-(1/(2*Cp))*(Dcj**2-gamma*RR*T2)
    eq3 = Dcj-(np.sqrt(gamma*RR*T0+qq*(gamma**2-1)/2)+np.sqrt(qq*(gamma**2-1)/2))
    #print(3*x[1]/x[0]) #effective activation energy
    #print(x)
    return [eq1,eq2,eq3]

def reflection(theta):#, P3,M_arr,d0,od):
    xx     = np.argmin(np.abs(r_l - (r_l[-1]-r0)))
    M_avg  = (M_l[-1])*math.sin(theta)#np.average(M_l[xx:-1])*math.sin(theta)
    P1     = ((1+(2*gamma)/(gamma+1)*(M_avg**2-1))*P0)
    rho1   = (((gamma+1)*M_avg**2)/(2+(gamma-1)*M_avg**2))*rho0
    T1     = (P1/(rho1*RR))
    R      = L/(2*math.cos(theta))
    def mach_ref(M_R):
        return (M_R/(M_R**2-1)) - (M_s/(M_s**2-1))*np.sqrt(1+(((2*(gamma-1))/(gamma+1)**2)*(M_s**2-1)*(gamma+1/M_s**2))) 
    xi     = np.argmin(np.abs(r_l - R))
    M_s    = M_l[xi]*math.sin(theta)*(np.sqrt(T0/T1))
    P2     = ((1+(2*gamma)/(gamma+1)*(M_s**2-1))*P1)
    rho2   = (((gamma+1)*M_s**2)/(2+(gamma-1)*M_s**2))*rho1
    T2     = (P2/(rho2*RR))
    M_R    = opt.root(mach_ref,M_s*0.95)
    PP3    = ((1+(2*gamma)/(gamma+1)*(M_R.x**2-1))*P2)
    rho3   = (((gamma+1)*M_s**2)/(2+(gamma-1)*M_s**2))*rho2
    #PP3    = q*rho3*((gamma-1)/(2*gamma-1))+pp3
    #print(R,'m_avg',M_avg, 'M_s = ',M_s,'M_R = ',M_R.x, T1,PP3)
    #print(M_s,rho2,P2,PP3,theta)
    return np.abs(PP3-P3)

def blast(AA,eff,T00,gam,qq,moll,P00,odd):
    ############################################
        #units SI
    global D0, C0, gamma,rho0,T0,P0,RR,Rr,alpha,h,A,Ea,q,rho_vn0, T_vn0, P_vn0, r_l, M_l, L, P3, od, r0, x0_star,r01,r02, t_ft
    vec       = [1500,2500,2000]
    x0_star   = 30
    r01       = 0.001
    r02       = 0.005
    t_ft      = 1.3
    #yvec      = [gam,moll,qq,T00]
    xx        = opt.fsolve(startup,vec,args=(gam,moll,qq,T00))
    A         = AA      # pre-exponential factor, m^3/kg-s
    Ea_hat    = eff*(xx[0]/T00)          # activation evergy, scaled by R*T0
    #print(Ea_hat)
    D0        = xx[2]         #m/s CJ velocity
    od        = odd    #1.785#1.7421#1.8764352818804846     #overdriven factor
    gamma     = gam          #specific heat ratio
    P0        = P00     #pa
    T0        = T00 
    mol       = moll        #  kg/mol
    Rr        = 8.3144087    # J/mol/K
    #tau0      = 1e-6        #characteristic reaction time
    #35.69444445  -330.62500005  1283.00694462 -2673.42708365  3160.74136143 -2016.17616683   545.06260004
    alpha     = 35.69444445*gamma**6-330.62500005*gamma**5+1283.00694462*gamma**4-2673.42708365*gamma**3+3160.74136143*gamma**2-2016.17616683*gamma+545.06260004   # point blast parameter  needs interpolation from table
    q         = qq*Rr*T0/mol           #currently not used because beta is chosen to cancel it
    h         = (np.pi)/4.0 # geometric factor 
    D_od      = D0*od #overdriven CJ velocity
    RR        = Rr/mol
    Ea        = Ea_hat*Rr*T0
    rho0      = P0/(RR*T0) #density
    C0        = np.sqrt(gamma*RR*T0)#sound speed in reactents mixture
    M0        = D_od/C0 
    #E0 = 50 # energy per length at point explision(d0) 
    beta      = 4*D0**2/(q*np.pi) #reaction parameter
    a_beta    = 4*D0**2
    rho_vn0    = (((gamma+1)*M0**2)/(2+(gamma-1)*M0**2))*rho0
    P_vn0      = ((1+(2*gamma)/(gamma+1)*(M0**2-1))*P0)
    T_vn0      = (P_vn0/(rho_vn0*RR))
    Mcj       = D0/C0
    P3        = (alpha*rho0*D_od**2*(gamma-1))/(h)
    #print('R0/d0 = ', np.sqrt(alpha*gamma)*Mcj*od,'\nMcj =',Mcj)
    #print(P3,P_vn0,alpha)
    ########################################
    #d0 =  opt.root(geo_mod,0.005)
    d0 =  opt.root_scalar(geo_mod,method='secant',x0=r01,x1=r02)
    #d0  =opt.bisect(geo_mod,0.0005,0.5)#,xtol=0.00001)
    #d0 = 0.005#0.005
    #r0        = d0.x/2
    r0        = d0.root/2
    ############################################
    t_star    = opt.fsolve(find_t_star,x0_star*r0/D0)#t_star    = ((np.sqrt(2)-1))*(od*r0/D0)opt.bisect(find_t_star,0.1*r0/D0,10*r0/D0)#
    r_star    = r(t_star)#(r0)*np.sqrt((1+od**2))
    t_end      = opt.fsolve(t_final,t_ft*t_star)

    time_jump  = t_star/1000
    print(t_end,t_star,t_end/t_star)
    t_arr      = np.arange(0, t_end+time_jump,time_jump) # t_arr      = np.arange(0,t_end+time_jump,time_jump)
    
    #############################################
    r_l       = np.fromiter(map(r, t_arr), dtype=np.float)
    D_l       = np.fromiter(map(D, t_arr), dtype=np.float)
    M_l       = np.fromiter(map(M, t_arr),dtype=np.float)
    #rho_l     = np.fromiter(map(rho_vn, t_arr),dtype=np.float)
    #P_l       = np.fromiter(map(P_vn, t_arr), dtype=np.float)
    #T_l       = np.fromiter(map(T_vn, t_arr), dtype=np.float)
    #tau_l     = np.fromiter(map(tau, t_arr), dtype=np.float)
    #u1_l      = np.fromiter(map(u1, t_arr), dtype=np.float)
    #x_ind_l   = np.fromiter(map(x_ind, t_arr), dtype=np.float)
    #x_ind2_l  = np.fromiter(map(x_ind2, t_arr), np.float)
    #x_ind3_l  = np.fromiter(map(x_ind3, t_arr), np.float)
    x_ind_r = x_ind_rel(t_star)
    #print(t_star,r_star,t_arr[-1],r_l[-1])
    ###################33

    #t_temp     = tau_l+t_arr
    #t_temp2    = tau_l+t_arrcell 
    # cell length is [cm] - 4.399900099981647 cell width is [cm] 1.8816314594251293 cell angle is [deg] 23.154129958901407
    #od is 1.5 r0 is [cm] - 0.13496455039336133

    #ix = np.abs(t_temp-t_arr[-1]).argmin()
    #t_temp =t_temp[:ix+1]
    #print(ix,np.abs((t_arr-t_star)).argmin(), t_star+tau_l[ix+1]-t_arr[-1],t_arr[ix]+tau_l[ix+1]-t_arr[-1])

    #r0_list.append(r0-(r_l[-1]-temp_front[ix]))

    #temp_front = (r_l) +(x_ind2_l)
    #temp_front = temp_front[:ix+1]

    #temp_front3 = (r_l) +(x_ind3_l)
    #temp_front3 = temp_front3[:ix+1]

    #temp_front2 = (r_l) +(x_ind_l)
    #temp_front2 = temp_front2[:ix+1]

    L = (r_l[-1]-r0)

    theta = opt.root(reflection, math.radians(30.0))
    L = L*100 #meters to centimeters
    W     = L*math.tan(theta.x)
    '''
    print('3 compression')
    print('!!!length units in centimeters!!!')
    print('cell length is [cm] -', L,' cell width is[cm] ', W,'theta is', math.degrees(theta.x))
    print('r0 is [cm] -', (r0)*100, 'new kernel size [cm]', (x_ind_r*100))
    print( 'D/D0 at last point' , D_l[-1]/D0,'\n')
    '''
    return[L,W,math.degrees(theta.x),od,r0*100]
'''
plt.plot(t_arr,r_l,label = 'pressure front')
plt.plot(t_temp,temp_front,label = 'temp front')
#plt.plot(t_arr[:ix+1],temp_front,label = 'temp front5')
plt.plot(t_temp,temp_front2,label = 'temp front2')
plt.plot(t_temp,temp_front3,label = 'temp front3')
plt.plot(t_arr,r0+D0*t_arr,label = 'const D')
plt.plot(t_star,r_star,'o',label = 'star')
plt.plot(t_arr[ix],r_l[ix],'o',label = 'last particle')
#plt.plot(t_temp,r_l[:ix+1])#x_ind2_l[:ix+1])#,labal='x_ind2')
#plt.plot(t_temp,x_ind2_l[:ix+1])#,labal='x_ind2')
#plt.plot(t_temp,ind_l[:ix+1])#,labal='x_ind2')
plt.grid()
plt.ylabel('r')
plt.xlabel('t')
plt.legend(loc='lower right',fontsize = 10)
plt.show()

D0_arr = np.full_like(t_arr,D0)
plt.plot(t_arr,D_l,label = 'pressure front')
#plt.plot(t,rtag,label = 'pressure front_tag')
plt.plot(t_arr,D0_arr,label = 'const D')
plt.grid()
plt.ylabel('D')
plt.xlabel('t')
plt.legend(loc='lower right',fontsize = 10)
plt.show()
'''
#print (r_arr,t_end)

