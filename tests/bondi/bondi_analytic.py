import numpy as np
from astropy import constants as const
from astropy import units as u
import matplotlib.pyplot as plt
import pdb


def get_Tfunc(T,r):
    utemp = C1 / (np.power(r, 2.) * np.power(T, n))
    result = (-2 / r + np.power(utemp, 2)) + (2. * (1. + n) * T + np.power((1. + n) * T, 2)) * (1. -2. / r + np.power(utemp, 2)) - C2prime
    return result

def define_globals(rs_in, mdot_in = 1., gam_in = 5./3):
    global rs, mdot, gam, C1, C2, C2prime
    rs = rs_in
    mdot = mdot_in
    gam = gam_in
    n = 1. / (gam - 1.)
    uc = np.sqrt(1 / (2. * rs))
    Vc = -np.sqrt(np.power(uc, 2.) / (1. - 3. * np.power(uc, 2)))
    Tc = -n * np.power(Vc, 2) / ((n + 1) * (n * np.power(Vc, 2) - 1.))
    C1 = uc * np.power(rs, 2) * np.power(Tc, n)
    C2 = np.power(1. + (1. + n) * Tc, 2.) * (1. - 2. / rs + np.power(C1, 2) / (np.power(rs, 4) * np.power(Tc, 2 * n)))
    uprime = C1 / (np.power(rs, 2) * np.power(Tc, n))
    C2prime = (-2. / rs + np.power(uprime, 2.)) + (2. * (1. + n) * Tc + np.power((1. + n) * Tc, 2))*(1. - 2. / rs + np.power(uprime, 2))

def get_T(r, ax = None, inflow_sol = True):
    rtol = 1.e-12
    ftol = 1.e-14
    Tinf = (np.sqrt(C2) - 1.) / (n + 1)
    Tapprox = np.power(C1 * np.sqrt(2. / np.power(r, 3)), 1. / n)
    
    bounds1 = [Tinf, Tapprox] # smaller T solution
    bounds2 = [np.fmax(Tapprox, Tinf), 1.] # larger T solution
    if inflow_sol:
        if (r < rs):
            Tmin = bounds1[0]
            Tmax = bounds1[1]
        else:
            Tmin = bounds2[0]
            Tmax = bounds2[1]
    else:
        if (r < rs):
            Tmin = bounds1[1]
            Tmax = 1.5
        else:
            Tmin = 1e-10 #0.1*Tinf
            Tmax = bounds2[0]

    if ax is not None:
        # this is just for test to visualize the function get_Tfunc()
        print("plotting")
        T_arr = np.linspace(Tmin, Tmax, 100)
        ax.plot(T_arr, get_Tfunc(T_arr, r),'k:')
        ax.axhline(0)
        ax.set_yscale('symlog')

    T0 = Tmin
    f0 = get_Tfunc(T0, r)
    T1 = Tmax
    f1 = get_Tfunc(T1, r)

    if (f0 * f1 > 0):
        print("error")
        return -1

    Th = 0.5*(T0 + T1)
    fh = get_Tfunc(Th, r)

    epsT = rtol*(Tmin + Tmax)


    while ((abs(Th - T0) > epsT) and (abs(Th - T1) > epsT)) and (abs(fh) > ftol):
        if (fh * f0 > 0): # bisection method
            T0 = Th
            f0 = fh
        else:
            T1 = Th
            f1 = fh

        Th = (T1 - T0) / 2. + T0
        fh = get_Tfunc(Th, r)
    
    return Th

def get_rho(T):
    Kn = 4 * np.pi * C1 / mdot
    n = 1. / (gam - 1.)
    return np.power(T, n) / Kn

def get_quantity_for_rarr(rarr, quantity, rs_in = np.power(10.,2.5), mdot_in = 1., gam_in = 5./3, inflow_sol=True):
    define_globals(rs_in, mdot_in, gam_in)
    Tarr = np.array([get_T(r, inflow_sol = inflow_sol) for r in rarr])
    rhoarr = get_rho(Tarr)
    if quantity =='T':
        return Tarr
    elif quantity == 'rho' or quantity == "RHO":
        return rhoarr
    elif quantity == 'ur' or quantity == 'u^r': 
        urarr = C1 / (np.power(rarr, 2) * np.power(Tarr, n))
        return urarr
    elif quantity == 'u' or quantity == "UU":
        uarr = rhoarr * Tarr * n
        return uarr
    else:
        return None

def _main():
    rarr= np.logspace(np.log10(2.1),np.log10(1e9),100)

    if 0:
        define_globals(1e5*3/16)
        #define_globals(np.sqrt(1e5))
        Tarr=[]
        for r in rarr:
            Tarr += [get_T(r,inflow_sol=True)] #C1,C2,n,

        Tarr= np.array(Tarr)
        rhoarr=np.power(Tarr,n)

        #plt.loglog(rarr,rhoarr) #*np.power(rarr,3./2))
        plt.loglog(rarr,Tarr*np.power(rarr,1./2))
        plt.savefig("./temp.png")
    else:
        # at one radius
        define_globals(16)#np.sqrt(1e5))
        #rb=80*256/(27*gam)
        T0=get_T(66.1135922573586)#100*rb)
        print(T0)
        #rho0 = np.power(T0,n)
        #r=50.980548842952516
        #rho = rho0*(r+rb)/r

        print("mdot = {:.18g}".format(4*np.pi*C1)) # to choose mdot
        fig,ax=plt.subplots(1,1,figsize=(8,6))
        T=get_T(1e8,ax,inflow_sol=True) #C1,C2,n,
        rho=np.power(T,n)
        u=rho*T*n
        print("T = {:.5g}, rho = {:.5g}, u = {:.5g}".format(T,rho,u))
        print("logrho = {:.6g}, loguoverrho = {:.6g}".format(np.log10(rho),np.log10(u/rho)))

if __name__ == "__main__":
    _main()
