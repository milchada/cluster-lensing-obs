#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from pylab import savefig, figure, show
from matplotlib.patches import Ellipse
from matplotlib.font_manager import FontProperties
import cosmolopy.distance as cd
import WLUtils.Cosmology as Cosmology

#===============================================#
#  Miyoung Choi                                 #
#  Mod: Oct. 2018                               #
#===============================================#


cosmo_par= Cosmology.setCosmology("WMAP7-ML")
cosmo = {'omega_M_0' : cosmo_par.Om0, 'omega_lambda_0' : cosmo_par.OL0, 'h' : cosmo_par.h}
cosmo = cd.set_omega_k_0(cosmo)
zd=0.2323
#zs=1.
Dd= cd.angular_diameter_distance(zd,0., **cosmo) #[Mpc]

print (np.pi*Dd/(10800.))

#gridsize = 0.14 * 60 / arcsecPerPix # number of pixels in a grid point ~ 170x170
gridsizeMpc = 14.4/60*(np.pi*Dd/(10800.)) # ~53Mpc
print gridsizeMpc

import glob, os
files = glob.glob(os.getcwd()+'/rhoproj/*cat.dat')
files.sort()

def main(filename):#=r'density_proj_2Gyr_cat.dat'
	snapnum = filename.split('proj_')[1].split('_cat')[0]

	x, y, z, g1, g2= np.genfromtxt(filename, unpack=True)

	fig= plt.figure()
	figtitle= 'Reduced Shear Field (without noise)'
	t = fig.text(0.5,0.95,figtitle,color='k',horizontalalignment='center',fontsize=14)
#===========direct plot of kappa map in the background=================#
	N=500
	ax=fig.add_subplot(111,aspect='equal')
	print x.min(), x.max(), y.min(), y.max()
	xi = np.linspace(x.min(), x.max(), N)
	yi = np.linspace(y.min(), y.max(), N)
	zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='linear')
	CK = ax.contourf(xi, yi, zi, N, cmap=plt.cm.Blues_r, alpha=None)
#=============generate shear field with same grid size as HST data binned===#

#xn=np.arange(x.min(), x.max(), gridsizeMpc)
#yn=np.arange(y.min(), y.max(), gridsizeMpc)
	xn=np.arange(x.min(), x.max(), 0.1)
	yn=np.arange(y.min(), y.max(), 0.1)
#xnArcmin=np.arange(x.min()/(np.pi*Dd/(10800.)), x.max()/(np.pi*Dd/(10800.)), pixsize)
#ynArcmin=np.arange(y.min()/(np.pi*Dd/(10800.)), y.max()/(np.pi*Dd/(10800.)), pixsize)
	G1= scipy.interpolate.griddata((x, y), g1, (xn[None,:], yn[:,None]), method='linear')
	G2= scipy.interpolate.griddata((x, y), g2, (xn[None,:], yn[:,None]), method='linear')
	K = scipy.interpolate.griddata((x, y), z, (xn[None,:], yn[:,None]), method='linear')
	for i in range(len(xn)):
  	    for j in range(len(xn)):
    		if K[i,j] < 0.1:
        	    gcplx=complex(G1[i,j],G2[i,j])/(1.-K[i,j])
        	    circles=[Ellipse((xn[j],yn[i]), width=50.*np.absolute(gcplx), \
                       height=0., angle=0.5*np.angle(gcplx,deg=True))]
            	    for e in circles:
          	    	ax.add_artist(e)
          		e.set_clip_box(ax.bbox)
          		e.set_facecolor('gold')
          		e.set_edgecolor('orange')
	plt.xlabel('[Mpc]')
	plt.ylabel('[Mpc]')
	plt.xlim(-2, 2)
	plt.ylim(-2, 2)
#plt.colorbar(CK, label='$ \Sigma [10^{10} M_{\odot} / Mpc^2] $', orientation='vertical')
	plt.colorbar(CK, label='$ \kappa $', orientation='vertical')
	plt.savefig('kappamaps/reducedshear_%s.png' % snapnum,dpi=300)
#plt.show()

for file in files:
	if files.index(file)>35:
		main(file)
