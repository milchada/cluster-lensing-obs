#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
from scipy import integrate
import Cosmology
import cosmolopy.distance as cd

#===============================================================================#
# Miyoung Choi                                                                  #
# Apr 2017                                                                      #
# Weak lensing stuff                                                            #
# Modifying/adding up to Brandyn's Cluster.py to handle CosmoOWLS fits files    #
#   - generate shear catalog from sigma maps                                    #
#       for sources at infinite redshift                                        #
#   - redshift distribution functions included                                  #
#   - class weakLens: beta calculator                                           #
#===============================================================================#


# Chang et al.
def pzfxn(z):
    alpha=1.24
    beta= 1.01
    z0=   0.51
    return (z**alpha)*np.exp(-(z/z0)**beta)
# integrated pzfxn for transformation method
def PZFT(z,fxn):
    I, err= integrate.quad(fxn,0.,z)
    return I

# general distribution function
def nzfxn(z):
    z0 = 2./3.
    return z**2./(2.*z0**3)*np.exp(-z/z0)

def NZFT(z,fxn):
    I, err= integrate.quad(fxn,0.,z)
    return I

# some astrophysical parameters
Mpc2meter=3.08568025*10**22 #Mpc to meters
clight=2.99792*10**8/Mpc2meter #m/s -> Mpc/s
#G = 6.67428*10**(-11)/(Mpc2meter**3) # m3/(kg·s^2) --> Mpc^3/(kg s^2)
G = 4.5177E-48 #units of Mpc^3/solar mass/s^2

class weakLens(object):
    def __init__(self, zd=None, cosmoName=None):

        if zd is None:
            print("No value set for the redshift. Setting z=0.25")
            self.zd = 0.25
        else:
            self.zd = zd

        if cosmoName is None:
            print("No cosmology provided. Using WMAP7-ML cosmology")
            cosmoName = "WMAP7-ML"
        self.cosmo = Cosmology.setCosmology(cosmoName)

    def beta(self, zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")
        cosmo = {'omega_M_0' : self.cosmo.Om0, \
                 'omega_lambda_0' : self.cosmo.OL0, \
                 'h' : self.cosmo.h}
        cosmo = cd.set_omega_k_0(cosmo)
        zd = self.zd
        Dd = cd.angular_diameter_distance(zd, **cosmo)
        Dcs = cd.comoving_distance(zs, z0=0., **cosmo)
        return 0. if zs<zd else 1.-Dd*(1.+zd)/Dcs

    def beta_inf(self):
        cosmo = {'omega_M_0' : self.cosmo.Om0, \
                 'omega_lambda_0' : self.cosmo.OL0, \
                 'h' : self.cosmo.h}
        cosmo = cd.set_omega_k_0(cosmo)
        zd = self.zd
        Dd = cd.angular_diameter_distance(zd, **cosmo)
        Dcs_inf = cd.comoving_distance(np.inf, z0=0., **cosmo)
        return 1.-Dd*(1.+zd)/Dcs_inf

    def sigmaC(self, zs=None):
        if zs is None:
            raise Exception("Error. Need to set beta")
        cosmo = {'omega_M_0' : self.cosmo.Om0, \
                 'omega_lambda_0' : self.cosmo.OL0, \
                 'h' : self.cosmo.h}
        cosmo = cd.set_omega_k_0(cosmo)
        zd = self.zd
        Dd = cd.angular_diameter_distance(zd, **cosmo)
        return clight**2/(4*np.pi*G*Dd*self.beta(zs))

