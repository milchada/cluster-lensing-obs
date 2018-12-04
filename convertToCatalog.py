import numpy as np
import sys, os
from WLUtils.weaklens import simLens as sl
from WLUtils.Cosmology import *
from astropy.io import fits
from astropy import units
import multiprocessing
#from joblib import Parallel, delayed

#===============================================#
#  Miyoung Choi                                 #
#  Mod: Oct. 2018                               #
#===============================================#


# set number of sigma maps
inputs= range(0,1)
num_cores=4

zd=0.2323
zs=1.

here=os.path.dirname(os.path.realpath(__file__))

def processInput(n, filename):
    #filename='density_proj_%dGyr.fits' % snapnum
    filename_out=filename.split('.fits')[0]+'_cat.dat'#density_proj_%dGyr_cat.dat'
    #if not os.path.exists(subdir):
    #  os.makedirs(subdir)
    if os.path.exists(os.path.join(here,filename)):
        hdulist=fits.open(filename)
        hdr = hdulist[0].header    
        print 'physical size per pixel = ', hdr[12], 'kpc'
        pixSize = hdr[12]
        
        sigma=fits.getdata(filename)
        conv = units.g.in_units('Msun')/(units.cm.in_units('kpc')**2)
	#print sigma, sigma.shape, sigma[0].shape, sigma[1].shape
        pixdim = len(sigma[0])
        fov = pixSize * float(len(sigma[0]))
        
        sigma=np.array(sigma)*conv*10**6 #convert the unit of density, Msun/kpc^2 --> Msun/Mpc^2
        
        lens = sl(sigma=sigma, z=zd, fovx=fov, fovy=fov)
        kappa          = sl.convergence(lens, zs=zs)
        gamma1, gamma2 = sl.shear(lens, zs=zs)



        #======save kappa and shear field======#
         
        xmpc = np.linspace(-fov/2./1000., fov/2./1000., pixdim)
        ympc = np.linspace(-fov/2./1000., fov/2./1000., pixdim)
        
        fw=open(os.path.join(here,filename_out),'w')
        for i in range(pixdim):
            for j in range(pixdim):
                fw.write('%0.10e %0.10e %0.10e %0.10e %0.10e' % (xmpc[j] \
                               , ympc[i] \
                               , kappa[i,j], gamma1[i,j], gamma2[i,j]))
                fw.write("\n")
        fw.close()
     
import glob
files = glob.glob('rhoproj/*fits')

pool = multiprocessing.Pool(num_cores)
for n in inputs:
	for filename in files:
		pool.apply_async(processInput, (n, filename))
#results = Parallel(n_jobs=num_cores)(delayed(processInput)(n) for n in inputs)
