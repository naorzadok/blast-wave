# znd 1 step (maybe a later version for multiple reaction will be made)
# All units are in SI

#main program
import numpy as np
import scipy.integrate as integ
import matplotlib as plt


def zndznd(gam,qq,aaa,Eaa,n,T00,P00,moll):
    #calculate timesteps, z vector [lam,x,u,rho]
    def proby(t,z):
        T = T0+(q/Cp)*z[0]+(1/(2*Cp))*(Dcj**2-z[2]**2)
        c      =np.sqrt(gamma*RR*T)
        dlam_dt=aa*((1-z[0])**n)*z[3]*np.exp(-Ea/(RR*T))
        #print(aa,z[0],Ea,RR,T)
        dx_dt  =z[2]  # x starts from 0
        du_dt  =((gamma-1)/(c**2-z[2]**2))*z[2]*q*dlam_dt #q scaled by RR*T0
        drho_dt=-(z[3]/z[2])*du_dt
        #dT_dt  =(q/Cp)*dlam_dt-u*du_dt 
        return [dlam_dt,dx_dt,du_dt,drho_dt]


    #file the program will write to
    file_name = 'test.txt'

    #set initaial conditions
    lam = 0.0            # mass fraction
    q   = qq     # heat of reaction, scaled by R*T0/mol
    aa  = aaa       # pre-exponential factor, m^3/kg-s
    Ea  = Eaa         # activation evergy, scaled by R*T0
    
    # heptane 1.29682198 17.14930528  0.02808893
    # C1 [ 1.2961726  17.18273912  0.02817085]

    #set problem conditions
    p0     = P00   # Pa
    T0     = T00    # K
    gamma  = gam
    mol    = moll     # kg/mol
    Rr     = 8.3144087 # J/mol/K
    RR     = Rr/mol

    f      = 1.0      # (D/Dcj)^2 ###
    x      = 0.0      # m (shock position)

    max_time = 5e-4
    max_timestep = 1.0e-7  # s
    #max_iteration = 500

    #  Calculate some inital state parameters...

    rho0   = p0/(RR*T0)
    e0     = p0/(rho0*(gamma-1.)) #probably internal energy, for enthalpy p0*gamma/(rho0*(gamma-1.))
    c0     = np.sqrt(gamma*RR*T0)

    q   = q  * RR * T0 
    Ea  = Ea * RR * T0

    #  Calculate ZND parameters as initial point for integration 

    qggg   = q/p0*rho0 * (gamma*gamma-1.) / (2.*gamma) #from gaznd
    Dcj    = c0 * (np.sqrt(1.+qggg) + np.sqrt(qggg) ) #from gaznd
    u0     = Dcj
    D      = np.sqrt(Dcj*Dcj*f)
    M1     = Dcj/c0 #mach number befora shockwave
    Cp     = (gamma*RR)/(gamma-1)
    p      = p0 * (1+((2*gamma)/(gamma+1)*(M1**2-1))) 
    rho    = rho0 *((M1**2*(gamma+1))/(2+(gamma-1)*M1**2)) 
    e      = e0 + 1/2. * (p+p0) * (1./rho0 - 1./rho)
    T1     = T0*(1+((2*(gamma-1))/(gamma+1)**2)*((gamma*M1**2+1)/(M1**2)*(M1**2-1)))
    u1     = u0*rho0/rho
    T_cj   = (2/(gamma+1))*(T0+q/Cp+(1/(2*Cp))*(Dcj**2))
    #p      = rho*RR*T1
    time   = 0.0
    c      = np.sqrt(gamma*RR*T1)
    T      = T1
    therm  = 0
    #breakout  = False #variable that turns True when close to sonic point
    
    eff_Ea = Eaa*T0/T1
    
    #write to file initial conditions
    #fo = open(file_name,'w')
    #fo.write('M1, p, rho , u1, T1, T_CJ= {:.5E} {:.5E} {:.5E} {:.5E} {:.5E} {:.5E} \n'.format(M1,p,rho,u1,T1,T_cj))
    #fo.write('ZND profile for Dcj, D = {:.5E} {:.5E} \n'.format(Dcj,D))
    #fo.write('VARIABLES= x,        p,           rho,          T,        U,          c,          M,         lam,          thermicity,         time \n') # lam1, w1, lam2, w2, time')

    #fo.write(x,p,rho,T1,u1,c,lam,therm,time)
    z = [lam,x,u1,rho]

    data = [[x],[p],[rho],[T],[u1],[c],[u1/c],[lam],[therm],[time]]
    #(x,p,rho,T1,u1,c,M,lam,therm,time) #array for other instances calculated 
    #after shock parameters,therm =0

    #sol = integ.solve_ivp(proby,t_span=(time,max_time),max_step=max_timestep,y0=z,method='DOP853',relTol=1e-5,absTol=1e-8,) #,args=(T) !t_eval=np.linspace(time,time+max_timestep,100


    sol = integ.solve_ivp(proby,t_span=(time,max_time),max_step=max_timestep,y0=z,method='RK45') #,args=(T) !t_eval=np.linspace(time,time+max_timestep,100
    #print(sol.y,sol.t)
    for j in range(0,len(sol.y[1])-1):  #calculate other parameters and save all parameters to data matrix
        #fo.write('  {:.6E}   {:.6E}    {:.6E}    {:.6E}    {:.6E}   {:.6E}   {:.6E}    {:.6E}    {:.6E}    {:.6E} \n'.format(data[0][-1],data[1][-1],data[2][-1],data[3][-1], data[4][-1],data[5][-1],data[6][-1],data[7][-1],data[8][-1],data[9][-1]))
        T = T0+(q/Cp)*sol.y[0][j+1]+(1/(2*Cp))*(Dcj**2-sol.y[2][j+1]**2)
        #print(sol.y[0][j+1])
        c = np.sqrt(gamma*RR*T)
        M = sol.y[2][j+1]/c
        if (c/sol.y[2][j+1]-1)<1e-3: #stop near the sonic point
            #fo.write('# sonic point arrived \n')
            #breakout = True
            break
        data[0].append(sol.y[1][j+1])
        data[1].append(sol.y[3][j+1]*RR*T)
        data[2].append(sol.y[3][j+1])
        data[3].append(T)
        data[4].append(sol.y[2][j+1])
        data[5].append(c)
        data[6].append(M)
        data[7].append(sol.y[0][j+1])
        data[8].append((sol.y[0][j+1]-sol.y[0][j])/(sol.t[j+1]-sol.t[j]))
        data[9].append(sol.t[j+1])
    #fo.write(data[0][-1],data[1][-1],data[2][-1],data[3][-1],data[4][-1],data[5][-1],data[6][-1],data[7][-1],data[8][-1],data[9][-1])
    #fo.write('# finished iterating')
    #fo.close()
    ind = data[8].index(max(data[8]))   #return index of max reaction rateind = data[8].index(max(data[8]))   #return index of max reaction rate
    Xdlist = [abs(i-0.5) for i in data[7]] 
    #print(Xdlist)
    Xd = Xdlist.index(min(Xdlist))
    #print('finished')
    return data, Xd ,eff_Ea
    #return Xd



