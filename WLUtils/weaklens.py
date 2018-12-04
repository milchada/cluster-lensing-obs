#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import Cosmology
import cosmolopy.distance as cd
from scipy import integrate



#===================================================================================#
# Miyoung Choi                                                                      #
# Latest Update: Apr 2017                                                           #
# Weak lensing stuff                                                                #
# Modifying/adding up to Brandyn's Cluster.py to handle suf.mass.den fits files     #
#   - fft to generate shear field                                                   #
#   - included beta                                                                 #
#   - rectangular field of view                                                     #
#===================================================================================#


# some astrophysical parameters
Mpc2meter=3.08568025*10**22 #Mpc to meters
clight=2.99792*10**8/Mpc2meter #m/s -> Mpc/s
#G = 6.67428*10**(-11)/(Mpc2meter**3) # m3/(kg·s^2) --> Mpc^3/(kg s^2)
G = 4.5177E-48 #units of Mpc^3/solar mass/s^2

#class for getting cosmo-OWLS cluster lens properties
class simLens(object):
    #initialize cluster. requires a surface density map, (solar masses/Mpc^2) redshift,
    # the field of view size (Mpc), and the cosmology in which it was generated.
    def __init__(self, sigma=None, z=None, fovx=None, fovy=None, cosmoName=None):
        if sigma is None:
            raise Exception("Error. Need to include 2D surface density map (solar masses/Mpc^2) as ndarray.")
        else:
            self.sigma = sigma
        if cosmoName is None:
            print("No cosmology provided. Using WMAP7-ML cosmology")
            cosmoName = "WMAP7-ML"
        self.cosmo = Cosmology.setCosmology(cosmoName)

        if z is None:
            print("No value set for the redshift. Setting z=0.25")
            self.z = 0.25
        else:
            self.z = z
            cosmo = {'omega_M_0' : self.cosmo.Om0, \
                 'omega_lambda_0' : self.cosmo.OL0, \
                 'h' : self.cosmo.h}
            cosmo = cd.set_omega_k_0(cosmo)
            self.Dd=cd.angular_diameter_distance(z, **cosmo) #Mpc

        if fovx is None:
            print("No field of view size set. Setting fov=30 Mpc")
            self.fovx = 30.
        else:
            self.fovx = fovx
            
        if fovy is None:
            print("No field of view size set. Setting fov=30 Mpc")
            self.fovy = 30.
        else:
            self.fovy = fovy

    # lensing depth 'beta' for sources at infinite redshift
    def beta_inf(self):
        cosmo = {'omega_M_0' : self.cosmo.Om0, \
                 'omega_lambda_0' : self.cosmo.OL0, \
                 'h' : self.cosmo.h}
        cosmo = cd.set_omega_k_0(cosmo)
        zd = self.z
        Dd = self.Dd
        Dcs_inf = cd.comoving_distance(np.inf, z0=0., **cosmo) #Mpc
        return 1.-Dd*(1.+zd)/Dcs_inf

    # lensing depth 'beta' to consider p(z). beta = 0 for z < zd
    def beta(self, zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")
        cosmo = {'omega_M_0' : self.cosmo.Om0, \
                 'omega_lambda_0' : self.cosmo.OL0, \
                 'h' : self.cosmo.h}
        cosmo = cd.set_omega_k_0(cosmo)
        zd = self.z
        Dd = self.Dd
        Dcs = cd.comoving_distance(zs, z0=0., **cosmo) #Mpc
        return 0. if zs<zd else 1.-Dd*(1.+zd)/Dcs


    #return critical surface density for cluster and source redshifts
    def sigmaC(self, zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")
        Dd = self.Dd
        return clight**2/(4*np.pi*G*Dd*self.beta(zs))

    def convergence(self,zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")
        return self.sigma/self.sigmaC(zs)

    def shear(self,zs=None, polar=False):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")

        #fft on convergence grid
        kappaF = np.fft.fft2(self.convergence(zs))
        pixx = self.fovx/kappaF.shape[0]
        pixy = self.fovy/kappaF.shape[1]

        #get the frequencies using pixel size as spacing
        freqx = np.fft.fftfreq(kappaF.shape[0],d=pixx)
        freqy = np.fft.fftfreq(kappaF.shape[1],d=pixy)
        freqX, freqY = np.meshgrid(freqx, freqy)
        #print kappaF.shape, freqX.shape
        #initialize and then calculate shear in fourier space
        gamma1F = np.zeros((kappaF.shape[0],kappaF.shape[1]), dtype=complex)
        gamma2F = np.zeros((kappaF.shape[0],kappaF.shape[1]), dtype=complex)

        gamma1F = kappaF*(freqX**2-freqY**2)/(freqX**2+freqY**2)
        gamma2F = 2.*kappaF*(freqX*freqY)/(freqX**2+freqY**2)

        #replace bad elements
        gamma1F = np.nan_to_num(gamma1F)
        gamma2F = np.nan_to_num(gamma2F)

        #fft back to position space
        gamma1 = np.fft.ifft2(gamma1F).real
        gamma2 = np.fft.ifft2(gamma2F).real

        return (gamma1, gamma2)

