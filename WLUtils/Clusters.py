import Cosmology,  math
import numpy as np
#from pymc import *
#
'''
class for clusters. Right now this only has spherical NFW and simulated clusters. Going to add triaxial later.

Use this for cluster mass properties. Now adding MCMC fitting.
'''

#some constants
G = 4.5177E-48 #units of Mpc^3/solar mass/s^2
clight = 9.716E-15 #units of Mpc/s

#use this class for analytical spherical NFW functions
class sphNFW(object):

    def __init__(self, c = None, r200 = None, M200 = None, z=None, cosmoName=None):

        if cosmoName is None:
            print("No cosmology provided. Using WMAP7-ML cosmology")
            cosmoName = "WMAP7-ML"
        self.cosmo = Cosmology.setCosmology(cosmoName)

        if z is None:
            print("No value set for the redshift. Setting z=0.25")
            self.z = 0.25
        else:
            self.z = z

        if c is None:
            print("No value set for concentration. Setting c = 4")
            self.c = 4.
        else:
            self.c = c

        if r200 is None and M200 is not None:
            self.M200 = M200
            self.r200 = self.M200Tor200()
        elif r200 is not None and M200 is None:
            self.r200 = r200
            self.M200 = self.r200ToM200()
        elif r200 is None and M200 is None:
            print ("No value set for r200 or M200. Setting r200=2 Mpc")
            self.r200 = 2.
            self.M200 = self.r200ToM200()
        else:
            raise Exception('You cannot set both M200 and r200')

    #convert r200 in Mpc to M200 in solar masses
    def r200ToM200(self):
        return 200*self.rhoC()*(4./3)*np.pi*(self.r200)**3

    #convert M200 in solar masses to r200 in Mpc
    def M200Tor200(self):
        return ((self.M200*3)/(200*self.rhoC()*4*np.pi))**(1./3)
    
    #return critical density at cluster redshift in units of solar masses/Mpc^3
    def rhoC(self):
        return self.cosmo.rho_c(self.z)*1000**3*self.cosmo.h**2

    #return delta_C for cluster concentration
    def deltaC(self):
        return (200./3.)*(self.c**3)/(np.log(1+self.c)-(self.c/(1+self.c)))

    #return critical surface density for cluster and source redshifts in solar masses/Mpc^2
    def sigmaC(self, zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")

        zl = self.z
        Dl = self.cosmo.angularDiameterDistance(zl) / self.cosmo.h
        Ds = self.cosmo.angularDiameterDistance(zs) / self.cosmo.h
        Dls = (1/(1+zs))*(Ds*(1+zs)-Dl*(1+zl))
        return clight**2*Ds/(4*math.pi*G*Dls*Dl)

    #return cluster scale radius in Mpc
    def rs(self):
        return self.r200 / self.c

    #Surface density profile. Input can be either 1D or 2D numpy ndarray or a float.
    def surfaceProfile(self, r=None):
        if r is None:
            raise Exception("Need to provide radius value or 1 or 2D numpy ndarray")

        rs = self.rs()
        k = 2*self.rhoC()*self.deltaC()

        if isinstance(r,np.ndarray):
            f = np.piecewise(r, [r>rs, r<rs,r==rs], [lambda r: (rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctan((((r/rs) - 1)/((r/rs) + 1))**0.5)/(((r/rs)**2) -1)**0.5)), lambda r:(rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctanh((((-r/rs) + 1)/((r/rs) + 1))**0.5)/((1 - ((r/rs)**2))**0.5))), lambda r: rs*k/3.])
        else:
            if r>rs:
                f = (rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctan((((r/rs) - 1)/((r/rs) + 1))**0.5)/(((r/rs)**2) -1)**0.5))
            elif r<rs:
                f = (rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctanh((((-r/rs) + 1)/((r/rs) + 1))**0.5)/((1 - ((r/rs)**2))**0.5)))
            else:
                f = rs*k/3.
        return f


    #calculate convergence on grid or as a profile. Function of radius in Mpc
    def convergence(self, r=None, zs=None):
        if r is None:
            raise Exception("Error. Need to provide radius value or 1 or 2D numpy ndarray")
        if zs is None:
            raise Exception("Error. Need to provide source redshift (zs)")

        return self.surfaceProfile(r) / self.sigmaC(zs)

    
    #calculate shear as a profile (1D numpy array) or single value. Function of radius in Mpc.
    def shear(self, r=None, zs=None):
        rs = self.rs()
        k = self.rhoC()*self.deltaC()
        sigmaC = self.sigmaC(zs)
        x = r/rs

        if isinstance(r,np.ndarray):
            f = np.piecewise(x, [x>1., x<1., x==1.], [lambda x: (rs*k)*((8*np.arctan(((x-1)/(1+x))**0.5)/(x**2*(x**2-1)**0.5)) \
                    + (4*np.log(x/2.)/x**2) \
                    - (2./(x**2-1)) \
                    + (4*np.arctan(((x-1)/(1+x))**0.5)/((x**2-1)**1.5))), \
                    lambda x: (rs*k)*((8*np.arctanh(((1-x)/(1+x))**0.5)/(x**2*(1-x**2)**0.5)) \
                    + (4*np.log(x/2.)/x**2) \
                    - (2./(x**2-1)) \
                    + (4*np.arctanh(((1-x)/(1+x))**0.5)/((x**2-1)*(1-x**2)**0.5))), \
                    lambda x: (rs*k)*(10./3+4.*np.log(0.5))])
        else:
            if (x<1.):
                f = (rs*k)*((8*np.arctanh(((1-x)/(1+x))**0.5)/(x**2*(1-x**2)**0.5)) \
                    + (4*np.log(x/2.)/x**2) \
                    - (2./(x**2-1)) \
                    + (4*np.arctanh(((1-x)/(1+x))**0.5)/((x**2-1)*(1-x**2)**0.5)))
            elif (x>1.):
                f = (rs*k)*((8*np.arctan(((x-1)/(1+x))**0.5)/(x**2*(x**2-1)**0.5)) \
                    + (4*np.log(x/2.)/x**2) \
                    - (2./(x**2-1)) \
                    + (4*np.arctan(((x-1)/(1+x))**0.5)/((x**2-1)**1.5)))
            else:
                f = (rs*k)*(10./3+4.*np.log(0.5))

        #print f/sigmaC
        return f/sigmaC

    #calculate reduced shear as a profile. Function of radius in Mpc. Returns tangential g
    def reducedShear(self,r=None,zs=None):
        gtan = self.shear(r,zs) / (1-self.convergence(r,zs))
        #print gmag
        return gtan

    #generate a catalog of galaxies at redshift zs with random locations and ellipticities lensed by spherical NFW cluster.
    #does not include any rejection techniques as of now. ndens is galaxy number density in arcmin^-2 and sige is dispersion of intrinsic ellip.
    #fov is the physical size of the grid over which you wanna make the catalog.
    def galCatGen(self, zs=None, ndens=None, sige=None, fov=None):
        if ndens is None:
            raise Exception("No galaxy number density given. Using 30 arcmin^-2")
        if sige is None:
            raise Exception("No intrinsic ellipticity dispersion given. Using 0.25")
        if fov is None:
            raise Exception("No fov given. Using 5 Mpc")

        #calculate the total number of galaxies in the fov
        dA = self.cosmo.angularDiameterDistance(zs) / (self.cosmo.h)
        radPerArcmin = np.pi/10800.
        MpcPerArcmin = dA*radPerArcmin #gives the number of Mpc per arcmin
        areaSqArcmin = fov**2/MpcPerArcmin**2 #area of fov in arcmin^2
        Ngal = ndens * areaSqArcmin #gives total number of galaxies in fov

        #generate random galaxy locations
        xphys = np.random.uniform(low = -fov/2.,high = fov/2., size = Ngal)
        yphys = np.random.uniform(low = -fov/2.,high = fov/2., size = Ngal)
        rphys = (xphys**2+yphys**2)**0.5

        gtan = self.reducedShear(r=rphys,zs=zs)
        epstan = gtan+np.random.normal(0.,sige,Ngal)

        return (rphys, epstan)


#class for getting cosmo-OWLS cluster lens properties
class simLens(object):
    #initialize cluster. requires a surface density map, (solar masses/Mpc^2) redshift,
    # the field of view size (Mpc), and the cosmology in which it was generated.
    def __init__(self, sigma=None, z=None, fov=None, cosmoName=None):
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

        if fov is None:
            print("No field of view size set. Setting fov=30 Mpc")
            self.fov = 30.
        else:
            self.fov = fov

    #return critical surface density for cluster and source redshifts in solar masses/Mpc^2
    def sigmaC(self, zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")

        zl = self.z
        Dl = self.cosmo.angularDiameterDistance(zl) / self.cosmo.h
        Ds = self.cosmo.angularDiameterDistance(zs) / self.cosmo.h
        Dls = (1/(1+zs))*(Ds*(1+zs)-Dl*(1+zl))
        return clight**2*Ds/(4*math.pi*G*Dls*Dl)

    #returns convergence map on grid
    def convergence(self,zs=None):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")
            
        return self.sigma/self.sigmaC(zs)

    #returns shear map grid calculated with FFT. Set polar=True to get (shear magnitude, shear orientation).
    #Otherwise returns (gamma1,gamma2).
    def shear(self,zs=None, polar=False):
        if zs is None:
            raise Exception("Error. Need to set source redshift (zs)")

        #fft on convergence grid
        kappaF = np.fft.fft2(self.convergence(zs))
        pix = self.fov/kappaF.shape[0]

        #get the frequencies using pixel size as spacing
        freqx = np.fft.fftfreq(kappaF.shape[0],d=pix)
        freqy = np.fft.fftfreq(kappaF.shape[1],d=pix)
        freqX, freqY = np.meshgrid(freqx, freqy)

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

        if polar:
            gammaMag = (gamma1**2+gamma2**2)**(1./2)
            gammaPhi = (1./2)*np.arctan2(gamma2,gamma1)
            return (gammaMag, gammaPhi)
        else:
            return (gamma1, gamma2)

    #calculate reduced shear on a grid for sim clusters. Set polar=True to get (g magnitude, g orientation).
    #Otherwise returns (g1,g2)
    def reducedShear(self,zs=None, polar=False, tanFlag=False):
        if zs is None:
            raise Exception("Error. Need to provide source redshift (zs)")

        gamma1,gamma2 = self.shear(zs)
        kappa = self.convergence(zs)
        g1 = gamma1/(1-kappa)
        g2 = gamma2/(1-kappa)

        if polar:
            gMag = (g1**2 + g2**2)**0.5
            gPhi = 0.5*np.arctan2(g2,g1)
            return (gMag,gPhi)
        else:
            return (g1,g2)

    #generate a catalog of galaxies at redshift zs with random locations and ellipticities lensed by cosmo-OWLS cluster.
    #does not include any rejection techniques as of now. ndens is galaxy number density in arcmin^-2 and sige is dispersion of intrinsic ellip.
    #set tanFlag to true to return only tangential component
    def galCatGen(self, zs=None, ndens=None, sige=None, tanFlag=False):
        if ndens is None:
            raise Exception("Error. No galaxy number density given.")
        if sige is None:
            raise Exception("Error. No intrinsic ellipticity dispersion given.")

        #calculate the total number of galaxies in the fov
        dA = self.cosmo.angularDiameterDistance(zs) / (self.cosmo.h)
        radPerArcmin = np.pi/10800.
        MpcPerArcmin = dA*radPerArcmin #gives the number of Mpc per arcmin
        areaSqArcmin = self.fov**2/MpcPerArcmin**2 #area of fov in arcmin^2
        Ngal = ndens * areaSqArcmin #gives total number of galaxies in fov

        #calculate reduced shear everywhere
        g1,g2 = self.reducedShear(zs=zs)
        pix = self.fov/g1.shape[0]
        npix = g1.shape[0]

        #generate random galaxy locations
        xind = np.random.randint(low = 0,high = npix, size = Ngal)
        yind = np.random.randint(low = 0, high = npix, size = Ngal)

        xphys = (xind*pix) - (self.fov/2.)
        yphys = (yind*pix) - (self.fov/2.)
        rphys = (xphys**2 + yphys**2)**0.5

        #get subset of g
        g1sub = g1[xind,yind]
        g2sub = g2[xind,yind]

        #add noise
        eps1 = g1sub+np.random.normal(0.,sige,Ngal)
        eps2 = g2sub+np.random.normal(0.,sige,Ngal)

        #print g1sub.shape
        #print xphys
        #print yphys

        if tanFlag:
            locPhi = np.arctan2(xphys,yphys)
            #print locPhi
            epstan = -(eps1*np.cos(2.*locPhi) + eps2*np.sin(2.*locPhi))
            return (rphys, epstan)
        else:
            return (rphys, eps1, eps2)
