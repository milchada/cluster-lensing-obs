from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import random
import multiprocessing
from joblib import Parallel, delayed

import os

#===============================================#
#  Miyoung Choi                                 #
#  Mod: Jun. 2018                               #
#===============================================#

#pixeldeg = 10./3600.
pixelmin = 10./60.
xarcmin = np.arange(-150., 150., pixelmin)
yarcmin = np.arange(-150., 150., pixelmin)
numgal = int(30*300*300)

def lookupNearest(x0, y0, data):
    xi = np.abs(xarcmin-x0).argmin()
    yi = np.abs(yarcmin-y0).argmin()
    data = data
    return data[yi,xi]

cosmology = 'WMAP9'
surveyName = 'LSST'
Mnus = [0.06, 0.12, 0.24, 0.48]
#Mnus = [0.06]
#Mnus = ['0_DMONLY', '0_low_AGN', '0_hi_AGN', '0']

inputs= range(0,25)

#for fi in range(25):
def processInput(fi):
#for fj in range(1,8):
  for Mnu in Mnus:
  #for AGN in AGNs:
    for i in range(5):
  
        #hdul = fits.open('WMAP9_CFHTLenS_rev_WL_maps/Mnu_0_DMONLY_WMAP9_cone_0_CFHTLenS_rev_3D_bin1.fits')
        #hdul = fits.open('WMAP9_KiDS_WL_maps/Mnu_0_'+AGN+'_WMAP9_cone_%.1d_KiDS_2D.fits'%fi)
        hdul = fits.open(cosmology+'_'+surveyName+'_2D_WL_maps/Mnu_'+str(Mnu)+'_'+cosmology+'_cone_%.1d_'%fi+surveyName+'_2D.fits')
        #hdul = fits.open('WMAP9_KiDS_WL_maps_no_star/Mnu_0_WMAP9_cone_%.1d_KiDS_2D_no_stars.fits'%fi)

        data = hdul[0].data # assuming/ the first extension is a table
        #print data.shape
        #print data[:,:,0].shape

        kappa = data[:,:,0]
        gamma1 = data[:,:,1]
        gamma2 = data[:,:,2]
        #print kappa.shape, gamma1.shape, gamma2.shape


        # dimensions: 5 deg x 5 deg with 10 arcsec pixels (1800x1800 pixels)




        '''
        #===============kappa map ============================#
            xv, yv = np.meshgrid(xarcmin,yarcmin)
        fig= plt.figure(0)
        ax=fig.add_subplot(111,aspect='equal')
        ax.contourf(xv, yv, kappa, np.linspace(kappa.min(),kappa.max(),100), cmap=plt.cm.binary_r, alpha=None)
        ax.set_xlabel('[deg]')
        ax.set_ylabel('[deg]')
        #plt.show()
        plt.savefig('BAHAMA_low_AGN_example_cone%.1d'%fi+'_bin%.1d.png'%fj, dpi=500)
        '''
        
        '''
        # full grid into cat
        f=open('WMAP9_KiDS_WL_cats/Mnu_0_low_AGN_WMAP9_cone_%.1df_KiDS_2D.txt'%fi, 'w')
        for i in range(1800):
            for j in range(1800):
                f.write('{0:10f} {1:10f} {2:10f} {3:10f}' \
                                # in arcmin unit
                                .format(xv[i,j]*60., \
                                        yv[i,j]*60., \
                                        gamma1[i,j], gamma2[i,j]))
                f.write("\n")
        f.close()
        '''
        

            
        # random point cat with the shape noise
        #f=open('WMAP9_KiDS_WL_cats/Mnu_0_'+AGN+'_WMAP9_cone_%.1d'%fi+'_KiDS_2D_realisation9_noise25_'+str(i)+'.txt', 'w')
        f=open(cosmology+'_'+surveyName+'_WL_cats/Mnu_'+str(Mnu)+'_'+cosmology+'_cone_%.1d'%fi+'_'+surveyName+'_2D_realisation30_noise25_'+str(i)+'.txt', 'w')
        #f=open('WMAP9_KiDS_WL_cats/Mnu_0_WMAP9_cone_%.1d'%fi+'_KiDS_2D_no_stars_realisation9_noise25_'+str(i)+'.txt', 'w')
        sig_eps=0.25
        en=np.random.normal(0.0,sig_eps,numgal)
        orient_s = np.pi*np.random.random(numgal)
        en1 = en*np.cos(2.*orient_s)
        en2 = en*np.sin(2.*orient_s)
        random.seed(fi)
        for igal in range(numgal):

            xnew = random.uniform(-150.,150.)
            ynew = random.uniform(-150.,150.)
            kappanew = (1.-lookupNearest(xnew, ynew, kappa))
            g1new = lookupNearest(xnew, ynew, gamma1)/(1.-lookupNearest(xnew, ynew, kappa))
            g2new = lookupNearest(xnew, ynew, gamma2)/(1.-lookupNearest(xnew, ynew, kappa))
            #print gam1new, gam2new
            gcplx = complex(g1new, g2new)
            escplx = complex(en1[igal], en2[igal])
            eicplx = (escplx+gcplx)/(1+gcplx.conjugate()*escplx)
            f.write('{0:10f} {1:10f} {2:10f} {3:10f}' \
                            # in arcmin unit
                            .format(xnew, ynew, eicplx.real, eicplx.imag))
            f.write("\n")
        f.close()
        print fi, 'th cat complete'


num_cores = multiprocessing.cpu_count()
pool = multiprocessing.Pool(num_cores)
for n in inputs:
	pool.apply_async(processInput, (n))
#results = Parallel(n_jobs=10)(delayed(processInput)(n) for n in inputs)
