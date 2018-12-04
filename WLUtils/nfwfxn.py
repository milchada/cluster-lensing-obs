#!`/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import cosmolopy.distance as cd
import Cosmology


def Sigma_nfw(x):
    #if isinstance(x,np.ndarray):
    if type(x) is np.ndarray:
        fnfw = np.piecewise(x, [x<1., x==1., x>1.], [lambda x:(1.-(2./np.sqrt(1.-x*x))*np.arctanh(np.sqrt((1.-x)/(1.+x))))/(x*x-1.), lambda x:1./3., lambda x:(1.-(2./np.sqrt(x*x-1.))*np.arctan2(np.sqrt(x-1.),np.sqrt(x+1.)))/(x*x-1.)])
        return fnfw

    else:
        if x<1:
            return (1.-(2./np.sqrt(1.-x*x))*np.arctanh(np.sqrt((1.-x)/(1.+x))))/(x*x-1.)
        elif x==1:
            return 1./3.
        else:
            return (1.-(2./np.sqrt(x*x-1.))*np.arctan2(np.sqrt(x-1.),np.sqrt(x+1.)))/(x*x-1.)
        


def Sigma_bar_nfw(x):
    #if isinstance(x,np.ndarray):

    if type(x) is np.ndarray:
        #print x
        fnfw = np.piecewise(x,[x<1., x==1., x>1.], \
                    [lambda x:(2./np.sqrt(1.-x*x)*np.arctanh(np.sqrt((1.-x)/(1.+x)))+np.log(x/2.))/(x*x), \
                     lambda x:(1.+np.log(0.5)), \
                     lambda x:(2./np.sqrt(x*x-1.)*np.arctan2(np.sqrt(x-1.),np.sqrt(1.+x))+np.log(x/2.))/(x*x)])
        return fnfw

    else:
        if x<1:
            return (2./np.sqrt(1.-x*x)*np.arctanh(np.sqrt((1.-x)/(1.+x)))+np.log(x/2.))/(x*x)
        elif x==1:
            return (1.+np.log(0.5))
        else:
            return (2./np.sqrt(x*x-1.)*np.arctan2(np.sqrt(x-1.),np.sqrt(1.+x))+np.log(x/2.))/(x*x)

    
# some astrophysical parameters
Mpc2meter=3.08568025*10**22 # Mpc to meters
Msolar2kg=1.989*10**30 # kg
clight=2.99792*10**8/Mpc2meter # m/s -> Mpc/s
G = 6.67428*10**(-11)/(Mpc2meter**3)*Msolar2kg # m3/(kg·s^2) --> Mpc^3/(solar mass * s^2)

# set cosmology
cosmo_par= Cosmology.setCosmology("WMAP7-ML")
cosmo = {'omega_M_0' : cosmo_par.Om0, \
         'omega_lambda_0' : cosmo_par.OL0, \
         'h' : cosmo_par.h}
cosmo = cd.set_omega_k_0(cosmo)

omegaM=cosmo_par.Om0
omegaL=cosmo_par.OL0
h=cosmo_par.h #no unit
H0=1000.*100.*h/Mpc2meter # /s
'''
# fac = (2.*rs*delta_c*rho_crit)/Sigma_crit
def facNFW(c, r200, zd, beta):
    #lens
    Dd= cd.angular_diameter_distance(zd,0., **cosmo)
    #print c, r200
    rs= r200/c
    delta_c=(200./3.)*(c**3)/(np.log(1.+c)-c/(1.+c))
    Sigma_crit = clight**2/(4.*np.pi*G*Dd*beta)

    rho_crit = 3.*(H0*H0)/(8.*np.pi*G)*(omegaM*(1.+zd)**3+omegaL) #solarmass/Mpc^3
    #print rs, delta_c, rho_crit, Sigma_crit
    #print (2.*rs*delta_c*rho_crit)/Sigma_crit
    return (2.*rs*delta_c*rho_crit)/Sigma_crit
'''
def kappa_nfw(r, c, r200, zd=None, beta=None):
    #def facNFW(c, r200, zd, beta):
    Dd= cd.angular_diameter_distance(zd,0., **cosmo)
    rs= r200/c
    delta_c=(200./3.)*(c**3)/(np.log(1.+c)-c/(1.+c))
    Sigma_crit = clight**2/(4.*np.pi*G*Dd*beta)
    rho_crit = 3.*(H0*H0)/(8.*np.pi*G)*(omegaM*(1.+zd)**3+omegaL) #solarmass/Mpc^3
    #return (2.*rs*delta_c*rho_crit)/Sigma_crit
    facNFW = (2.*rs*delta_c*rho_crit)/Sigma_crit
    return facNFW * Sigma_nfw(r/rs)

def gamma_nfw(r, c, r200, zd=None, beta=None): #kappa_bar - kappa
    #def facNFW(c, r200, zd, beta):
    #lens
    Dd= cd.angular_diameter_distance(zd,0., **cosmo)
    #print c, r200
    rs= r200/c
    delta_c=(200./3.)*(c**3)/(np.log(1.+c)-c/(1.+c))
    Sigma_crit = clight**2/(4.*np.pi*G*Dd*beta)
    rho_crit = 3.*(H0*H0)/(8.*np.pi*G)*(omegaM*(1.+zd)**3+omegaL) #solarmass/Mpc^3
    #return (2.*rs*delta_c*rho_crit)/Sigma_crit
    facNFW = (2.*rs*delta_c*rho_crit)/Sigma_crit
    return 2.*facNFW * Sigma_bar_nfw(r/rs) - facNFW * Sigma_nfw(r/rs)

def g_nfw(r, c, r200, zd=None, beta=None):
    # r[Mpc], c, r200, zd, beta
    Dd= cd.angular_diameter_distance(zd,0., **cosmo)
    #print c, r200
    rs= r200/c
    delta_c=(200./3.)*(c**3)/(np.log(1.+c)-c/(1.+c))
    Sigma_crit = clight**2/(4.*np.pi*G*Dd*beta)
    rho_crit = 3.*(H0*H0)/(8.*np.pi*G)*(omegaM*(1.+zd)**3+omegaL) #solarmass/Mpc^3
    #return (2.*rs*delta_c*rho_crit)/Sigma_crit
    facNFW = (2.*rs*delta_c*rho_crit)/Sigma_crit
    #return gamma_nfw(x, c, r200, zd, beta)/(1.-kappa_nfw(x, c, r200, zd, beta))
    g = (2.*facNFW * Sigma_bar_nfw(r/rs) - facNFW * Sigma_nfw(r/rs))/(1.-facNFW * Sigma_nfw(r/rs))
    #print type(x)
    return g

def mue(x, c, r200, zd=None, beta=None):
    return 1./np.absolute((1.-kappa_nfw(x, c, r200, zd, beta))**2-(gamma_nfw(x, c, r200, zd, beta)**2))

def mueinvsqrt(r, c, r200, zd=None, beta=None):
    #return np.sqrt(np.absolute((1.-kappa_nfw(x,fac))**2-(gamma_nfw(x,fac)**2)))
    Dd= cd.angular_diameter_distance(zd,0., **cosmo)
    rs= r200/c
    delta_c=(200./3.)*(c**3)/(np.log(1.+c)-c/(1.+c))
    Sigma_crit = clight**2/(4.*np.pi*G*Dd*beta)
    rho_crit = 3.*(H0*H0)/(8.*np.pi*G)*(omegaM*(1.+zd)**3+omegaL) #solarmass/Mpc^3
    #return (2.*rs*delta_c*rho_crit)/Sigma_crit
    facNFW = (2.*rs*delta_c*rho_crit)/Sigma_crit
    return np.sqrt(1./(1.+2.*facNFW * Sigma_nfw(r/rs)))

def mueIntegrand(x, c, r200, zd=None, beta=None):
    return x*np.sqrt(np.absolute((1.-kappa_nfw(x, c, r200, zd, beta))**2-(gamma_nfw(x, c, r200, zd, beta)**2)))

