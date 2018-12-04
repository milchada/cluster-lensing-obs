import numpy as np
import Cosmology
from weaklens import simLens as sl
from astropy.io import fits
import multiprocessing
from joblib import Parallel, delayed
import os, sys, csv, operator


inputs= range(0,1000)
#num_cores = multiprocessing.cpu_count()
num_cores=4

zd=0.25
zs=1.
pixdim=2000
here=os.path.dirname(os.path.realpath(__file__))

def processInput(n):
  filename="cluster_%.1d_z.fits"%n  
  filename_out="cluster_zd25_zs1_%.1d_3x3Mpc2arcmin.dat"%n
  subdir="clusters_zd25_zs1_3x3Mpc2arcmin"
  if not os.path.exists(subdir):
    os.makedirs(subdir)
  if os.path.exists(os.path.join(here,filename)):
    hdulist=fits.open(filename)
    sigma=fits.getdata(filename)
    sigma=np.array(sigma)
    sample=sl(sigma=sigma, z=zd)
    kappa          = sl.convergence(sample,zs=zs)
    gamma1, gamma2 = sl.shear(sample,zs=zs)
    print kappa.shape, gamma1.shape
    fw=open(os.path.join(here,subdir,filename_out),'w')
  
    #xmpc=np.linspace(-15.,15.,pixdim)
    #ympc=np.linspace(-15.,15.,pixdim)
    Dd=sample.Dd
    print n, Dd
    xarcmin=np.linspace(-15./(np.pi*Dd/(10800.)),15./(np.pi*Dd/(10800.)),pixdim)
    yarcmin=np.linspace(-15./(np.pi*Dd/(10800.)),15./(np.pi*Dd/(10800.)),pixdim)
    for i in range(pixdim):
      for j in range(pixdim):
      #print i*pixdim+j
      # i row (y-dir), j column (x-dir)
        if xarcmin[i] <= 3./(np.pi*Dd/(10800.)) \
           and xarcmin[i] >= -3./(np.pi*Dd/(10800.)) \
           and yarcmin[j] <= 3./(np.pi*Dd/(10800.)) \
           and yarcmin[j] >= -3./(np.pi*Dd/(10800.)):
      # writing position in unit of arcmin
           fw.write('{0:10f} {1:10f} {2:10f} {3:10f} {4:10f}' \
                 .format(xarcmin[i] \
                       , yarcmin[j] \
                       , kappa[j,i], gamma1[i,j], gamma2[i,j])
                )
           fw.write("\n")
    fw.close()
 
  

results = Parallel(n_jobs=num_cores)(delayed(processInput)(n) for n in inputs)
