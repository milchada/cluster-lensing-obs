###################################################################################################
#
# Cosmology.py 		(c) Benedikt Diemer
#						University of Chicago
#     				    bdiemer@oddjob.uchicago.edu
#
###################################################################################################

"""
This module is an implementation of the standard :math:`\Lambda CDM` cosmology, with a focus on 
structure formation applications. It assumes a fixed dark energy equation of state, w = -1, and 
includes the contributions from dark matter, dark energy, baryons, curvature, photons, and
neutrinos.

---------------------------------------------------------------------------------------------------
Basic usage
---------------------------------------------------------------------------------------------------

***************************************************************************************************
Setting and getting cosmologies
***************************************************************************************************

Colossus does not set a default cosmology. Thus, before using any cosmological functions, or any 
other Colossus functions that use the Cosmology module, a cosmology needs to be set. This is almost 
always done with the :func:`setCosmology` function, which can be used in multiple ways:

* Set one of the pre-defined cosmologies::
	
	setCosmology('WMAP9')

* Set one of the pre-defined cosmologies, but overwrite certain parameters::
	
	setCosmology('WMAP9', {'print_warnings': False})

* Add a new cosmology to the global list of available cosmologies. This has the advantage that the 
  new cosmology can be set from anywhere in the code. Only the main cosmological parameters are 
  mandatory, all other parameters can be left to their default values::
	
	params = {'flat': True, 'H0': 67.2, 'Om0': 0.31, 'Ob0': 0.049, 'sigma8': 0.81, 'ns': 0.95}
	addCosmology('myCosmo', params)
	cosmo = setCosmology('myCosmo')

* Set a new cosmology without adding it to the global list of available cosmologies::
	
	params = {'flat': True, 'H0': 67.2, 'Om0': 0.31, 'Ob0': 0.049, 'sigma8': 0.81, 'ns': 0.95}
	cosmo = setCosmology('myCosmo', params)

* Set a self-similar cosmology with a power-law power spectrum of a certain slope, and the 
  default settings set in the ``powerlaw`` cosmology::
	
	cosmo = setCosmology('powerlaw_-2.60')

Whichever way a cosmology is set, the current cosmology is stored in a global variable and 
can be obtained at any time::
	
	cosmo = getCurrent()

***************************************************************************************************
Changing and switching cosmologies
***************************************************************************************************

The current cosmology can also be set to an already existing cosmology object, for example when
switching between cosmologies::

	cosmo1 = setCosmology('WMAP9')
	cosmo2 = setCosmology('planck13')
	setCurrent(cosmo1)

The user can change the cosmological parameters of an existing cosmology object at run-time, but 
MUST call the update function directly after the changes. This function ensures that the parameters 
are consistent (e.g., flatness), and discards pre-computed quantities::

	cosmo = setCosmology('WMAP9')
	cosmo.Om0 = 0.31
	cosmo.checkForChangedCosmology()

***************************************************************************************************
Summary of getter and setter functions
***************************************************************************************************

.. autosummary::
	setCosmology
	addCosmology
	setCurrent
	getCurrent

***************************************************************************************************
Standard cosmologies
***************************************************************************************************

============== ===================== =========== =======================================
ID             Paper                 Location    Explanation
============== ===================== =========== =======================================
planck15-only  Planck Collab. 2015   Table 4     Best-fit, Planck only (column 2) 					
planck15       Planck Collab. 2015 	 Table 4     Best-fit with ext (column 6)			
planck13-only  Planck Collab. 2013   Table 2     Best-fit, Planck only 					
planck13       Planck Collab. 2013 	 Table 5     Best-fit with BAO etc. 					
WMAP9-only     Hinshaw et al. 2013   Table 1     Max. likelihood, WMAP only 				
WMAP9-ML       Hinshaw et al. 2013   Table 1     Max. likelihood, with eCMB, BAO and H0 	
WMAP9          Hinshaw et al. 2013   Table 4     Best-fit, with eCMB, BAO and H0 		
WMAP7-only     Komatsu et al. 2011   Table 1     Max. likelihood, WMAP only 				
WMAP7-ML       Komatsu et al. 2011   Table 1     Max. likelihood, with BAO and H0 		
WMAP7 	       Komatsu et al. 2011   Table 1     Best-fit, with BAO and H0 				
WMAP5-only     Komatsu et al. 2009   Table 1     Max. likelihood, WMAP only 			
WMAP5-ML       Komatsu et al. 2009   Table 1     Max. likelihood, with BAO and SN 		
WMAP5 	       Komatsu et al. 2009   Table 1     Best-fit, with BAO and SN 			
WMAP3-ML       Spergel et al. 2007   Table 2     Max.likelihood, WMAP only 				
WMAP3          Spergel et al. 2007   Table 5     Best fit, WMAP only 					
WMAP1-ML       Spergel et al. 2005   Table 1/4   Max.likelihood, WMAP only 				
WMAP1          Spergel et al. 2005   Table 7/4   Best fit, WMAP only 					
bolshoi	       Klypin et al. 2011    --          Cosmology of the Bolshoi simulation
millennium     Springel et al. 2005	 --          Cosmology of the Millennium simulation 
powerlaw       --                    --          Default settings for power-law cosms.
============== ===================== =========== =======================================

Those cosmologies that refer to particular simulations (such as bolshoi and millennium) are set 
to ignore relativistic species, i.e. photons and neutrinos.

***************************************************************************************************
Units
***************************************************************************************************

Unless otherwise noted, all functions in this module use the following units:

============ =======================================
Variable     Unit
============ =======================================
Length       Comoving Mpc/h
Wavenumber   Comoving h/Mpc
Time         Gigayears
Mass         :math:`M_{\odot}/h`
Density      Physical :math:`M_{\odot} h^2 / kpc^3`
============ =======================================

Note that densities use physical kpc/h whereas all other lengths are expressed in comoving Mpc/h. 
The kpc/h units for density make them compatible with the HaloDensityProfile module, where all 
lengths are expressed in physical kpc/h.

***************************************************************************************************
Derivatives and Inverses
***************************************************************************************************

Almost all cosmology functions that are interpolated (e.g., :func:`Cosmology.age()`, 
:func:`Cosmology.luminosityDistance()` or :func:`Cosmology.sigma()`) can be evaluated as an nth 
derivative. Please note that some functions are interpolated in log space, resulting in a logarithmic
derivative, while others are interpolated and differentiated in linear space. Please see the 
function documentations below for details.

The derivative functions were not systematically tested for accuracy. Their accuracy will depend
on how well the function in question is represented by the spline approximation. In general, 
the accuracy of the derivatives will be worse that the error quoted on the function itself, and 
get worse with the order of the derivative.

Furthermore, the inverse of interpolated functions can be evaluated by passing ``inverse = True``.
In this case, for a function y(x), x(y) is returned instead. Those functions raise an Exception if
the requested value lies outside the range of the interpolating spline.

The inverse and derivative flags can be combined to give the derivative of the inverse, i.e. dx/dy. 
Once again, please check the function documentation whether that derivative is in linear or 
logarithmic units.

---------------------------------------------------------------------------------------------------
Performance optimization and accuracy
---------------------------------------------------------------------------------------------------

This module is optimized for fast performance, particularly in computationally intensive
functions such as the correlation function. Almost all quantities are, by 
default, tabulated, stored in files, and re-loaded when the same cosmology is set again. For 
some rare applications (for example, MCMC chains where functions are evaluated few times, but for 
a large number of cosmologies), the user can turn this behavior off::

	cosmo = Cosmology.setCosmology('WMAP9', {"interpolation": False, "storage": False})

For more details, please see the documentation of the ``interpolation`` and ``storage`` parameters.
In order to turn off the interpolation temporarily, the user can simply switch the ``interpolation``
parameter off::
	
	cosmo.interpolation = False
	Pk = cosmo.matterPowerSpectrum(k)
	cosmo.interpolation = True
	
In this example, the power spectrum is evaluated directly without interpolation. The 
interpolation is accurate to better than 0.2% unless specifically noted in the function 
documentation, meaning that it is very rarely necessary to use the exact routines. 

---------------------------------------------------------------------------------------------------
Detailed Documentation
---------------------------------------------------------------------------------------------------
"""

###################################################################################################

import os
import math
import numpy
import scipy.integrate
import scipy.special
import scipy.interpolate
import hashlib
import pickle
import Utilities

###################################################################################################
# Useful constants
###################################################################################################

AST_pc  = 3.08568025E18
"""A parsec in centimeters."""
AST_kpc = 3.08568025E21
"""A kiloparsec in centimeters."""
AST_Mpc = 3.08568025E24 
"""A megaparsec in centimeters."""
AST_year = 31556926.0
"""A year in seconds."""
AST_Msun = 1.98892E33
"""A solar mass, :math:`M_{\odot}`, in grams."""
AST_c = 2.99792458E10
"""The speed of light in cm/s."""
AST_G = 4.30172E-6
"""The gravitational constant G in :math:`kpc \ km^2 / M_{\odot} / s^2`."""
AST_rho_crit_0_kpc3 = 2.774848e+02
"""The critical density of the universe at z = 0 in units of :math:`M_{\odot} h^2 / kpc^3`."""
AST_rho_crit_0_Mpc3 = 2.774848e+11
"""The critical density of the universe at z = 0 in units of :math:`M_{\odot} h^2 / Mpc^3`."""
AST_delta_collapse = 1.686
"""The threshold overdensity for halo collapse according to the top-hat collapse model."""

###################################################################################################
# Global variables for cosmology object and pre-set cosmologies
###################################################################################################

# This variable should never be used by the user directly, but instead be handled with getCurrent
# and setCosmology.
current_cosmo = None

# The following named cosmologies can be set by calling setCosmology(name). Note that changes in
# cosmological parameters are tracked to the fourth digit, which is why all parameters are rounded
# to at most four digits. See documentation at the top of this file for references.
cosmologies = {}
cosmologies['planck15-only'] = {'flat': True, 'H0': 67.81, 'Om0': 0.3080, 'Ob0': 0.0484, 'sigma8': 0.8149, 'ns': 0.9677}
cosmologies['planck15']      = {'flat': True, 'H0': 67.74, 'Om0': 0.3089, 'Ob0': 0.0486, 'sigma8': 0.8159, 'ns': 0.9667}
cosmologies['planck13-only'] = {'flat': True, 'H0': 67.11, 'Om0': 0.3175, 'Ob0': 0.0490, 'sigma8': 0.8344, 'ns': 0.9624}
cosmologies['planck13']      = {'flat': True, 'H0': 67.77, 'Om0': 0.3071, 'Ob0': 0.0483, 'sigma8': 0.8288, 'ns': 0.9611}
cosmologies['WMAP9-only']    = {'flat': True, 'H0': 69.70, 'Om0': 0.2814, 'Ob0': 0.0464, 'sigma8': 0.8200, 'ns': 0.9710}
cosmologies['WMAP9-ML']      = {'flat': True, 'H0': 69.70, 'Om0': 0.2821, 'Ob0': 0.0461, 'sigma8': 0.8170, 'ns': 0.9646}
cosmologies['WMAP9']         = {'flat': True, 'H0': 69.32, 'Om0': 0.2865, 'Ob0': 0.0463, 'sigma8': 0.8200, 'ns': 0.9608}
cosmologies['WMAP7-only']    = {'flat': True, 'H0': 70.30, 'Om0': 0.2711, 'Ob0': 0.0451, 'sigma8': 0.8090, 'ns': 0.9660}
cosmologies['WMAP7-ML']      = {'flat': True, 'H0': 70.40, 'Om0': 0.2715, 'Ob0': 0.0455, 'sigma8': 0.8100, 'ns': 0.9670}
cosmologies['WMAP7']         = {'flat': True, 'H0': 70.20, 'Om0': 0.2743, 'Ob0': 0.0458, 'sigma8': 0.8160, 'ns': 0.9680}
cosmologies['WMAP5-only']    = {'flat': True, 'H0': 72.40, 'Om0': 0.2495, 'Ob0': 0.0432, 'sigma8': 0.7870, 'ns': 0.9610}
cosmologies['WMAP5-ML']      = {'flat': True, 'H0': 70.20, 'Om0': 0.2769, 'Ob0': 0.0459, 'sigma8': 0.8170, 'ns': 0.9620}
cosmologies['WMAP5']         = {'flat': True, 'H0': 70.50, 'Om0': 0.2732, 'Ob0': 0.0456, 'sigma8': 0.8120, 'ns': 0.9600}
cosmologies['WMAP3-ML']      = {'flat': True, 'H0': 73.20, 'Om0': 0.2370, 'Ob0': 0.0414, 'sigma8': 0.7560, 'ns': 0.9540}
cosmologies['WMAP3']         = {'flat': True, 'H0': 73.50, 'Om0': 0.2342, 'Ob0': 0.0413, 'sigma8': 0.7420, 'ns': 0.9510}
cosmologies['WMAP1-ML']      = {'flat': True, 'H0': 68.00, 'Om0': 0.3136, 'Ob0': 0.0497, 'sigma8': 0.9000, 'ns': 0.9700}
cosmologies['WMAP1']         = {'flat': True, 'H0': 72.00, 'Om0': 0.2700, 'Ob0': 0.0463, 'sigma8': 0.9000, 'ns': 0.9900}
cosmologies['bolshoi']       = {'flat': True, 'H0': 70.00, 'Om0': 0.2700, 'Ob0': 0.0469, 'sigma8': 0.8200, 'ns': 0.9500, 'relspecies': False}
cosmologies['millennium']    = {'flat': True, 'H0': 73.00, 'Om0': 0.2500, 'Ob0': 0.0450, 'sigma8': 0.9000, 'ns': 1.0000, 'relspecies': False}
cosmologies['powerlaw']      = {'flat': True, 'H0': 70.00, 'Om0': 1.0000, 'Ob0': 0.0000, 'sigma8': 0.8200, 'ns': 1.0000, 'relspecies': False}

###################################################################################################
# Cosmology class
###################################################################################################

class Cosmology(object):
	"""
	A cosmology is set through the parameters passed to the constructor. Any parameter whose default
	value is ``None`` must be set by the user. This can easily be done using the 
	:func:`setCosmology()` function with one of the pre-defined sets of cosmological parameters 
	listed above. 
	
	Some parameters that are well constrained and have a sub-dominant impact on the computations
	have pre-set default values, such as the CMB temperature (T = 2.7255 K) and the effective number 
	of neutrino species (Neff = 3.046). These values are compatible with the most recent 
	measurements and can be changed by the user. 
	
	Parameters
	-----------------------------------------------------------------------------------------------
	name: str		
		A name for the cosmology, e.g. ``WMAP9``.
	flat: bool
		If flat, there is no curvature, :math:`\Omega_k = 0`, and :math:`\Omega_{\Lambda} = 1 - \Omega_m`.
	relspecies: bool
		If False, all relativistic contributions to the energy density of the universe (such as 
		photons and neutrinos) are ignored.
	Om0: float
		:math:`\Omega_m` at z = 0.
	OL0: float
		:math:`\Omega_{\Lambda}` at z = 0. This parameter is ignored if ``flat == True``.
	Ob0: float
		:math:`\Omega_{baryon}` at z = 0.
	H0: float
		The Hubble constant in km/s/Mpc.
	sigma8: float
		The normalization of the power spectrum, i.e. the variance when the field is filtered with a 
		top hat filter of radius 8 Mpc/h.
	ns: float
		The tilt of the primordial power spectrum.
	Tcmb0: float
		The temperature of the CMB today in Kelvin.
	Neff: float
		The effective number of neutrino species.
	power_law: bool
		Assume a power-law matter power spectrum, :math:`P(k) = k^{power\_law\_n}`.
	power_law_n: float
		See ``power_law``.
	interpolation: bool
		By default, lookup tables are created for certain computationally intensive quantities, 
		cutting down the computation times for future calculations. If ``interpolation == False``,
		all interpolation is switched off. This can be useful when evaluating quantities for many
		different cosmologies (where computing the tables takes a prohibitively long time). 
		However, many functions will be *much* slower if this setting is False, please use it only 
		if absolutely necessary. Furthermore, the derivative functions of :math:`P(k)`, 
		:math:`\sigma(R)` etc will not work if ``interpolation == False``.
	storage: bool 
		By default, interpolation tables and other data are stored in a permanent file for
		each cosmology. This avoids re-computing the tables when the same cosmology is set again. 
		However, if any file access is to be avoided (for example in MCMC chains), the user can 
		set ``storage = False``.
	print_info: bool
		Output information to the console.
	print_warnings: bool
		Output warnings to the console.
	text_output:
		If True, all persistent data (such as lookup tables for :math:`\sigma(R)` etc) is written 
		into named text files in addition to the default storage system. This feature allows the 
		use of these tables outside of this module. *Warning*: Be careful with the text_output 
		feature. Changes in cosmology are not necessarily reflected in the text file names, and 
		they can thus overwrite the correct values. Always remove the text files from the directory 
		after use. 
	"""
	
	def __init__(self, name = None, \
		Om0 = None, OL0 = None, Ob0 = None, H0 = None, sigma8 = None, ns = None, \
		flat = True, relspecies = True, Tcmb0 = 2.7255, Neff = 3.046, \
		power_law = False, power_law_n = 0.0, \
		print_info = False, print_warnings = True, \
		interpolation = True, storage = True, text_output = False):
		
		if name is None:
			raise Exception('A name for the cosmology must be set.')
		if Om0 is None:
			raise Exception('Parameter Om0 must be set.')
		if Ob0 is None:
			raise Exception('Parameter Ob0 must be set.')
		if H0 is None:
			raise Exception('Parameter H0 must be set.')
		if sigma8 is None:
			raise Exception('Parameter sigma8 must be set.')
		if ns is None:
			raise Exception('Parameter ns must be set.')
		if Tcmb0 is None:
			raise Exception('Parameter Tcmb0 must be set.')
		if Neff is None:
			raise Exception('Parameter Neff must be set.')
		if power_law and power_law_n is None:
			raise Exception('For a power-law cosmology, power_law_n must be set.')
		if not flat and OL0 is None:
			raise Exception('OL0 must be set for non-flat cosmologies.')
	
		# Copy the cosmological parameters into the class
		self.name = name
		self.flat = flat
		self.relspecies = relspecies
		self.power_law = power_law
		self.power_law_n = power_law_n
		self.Om0 = Om0
		self.OL0 = OL0
		self.Ob0 = Ob0
		self.H0 = H0
		self.sigma8 = sigma8
		self.ns = ns
		self.Tcmb0 = Tcmb0
		self.Neff = Neff

		# Compute some derived cosmological variables
		self.h = H0 / 100.0
		self.h2 = self.h**2
		self.Omh2 = self.Om0 * self.h2
		self.Ombh2 = self.Ob0 * self.h2
		
		if self.relspecies:
			# To convert the CMB temperature into a fractional energy density, we follow these
			# steps:
			# 
			# rho_gamma   = 4 sigma_SB / c * T_CMB^4 [erg/cm^3]
			#             = 4 sigma_SB / c^3 * T_CMB^4 [g/cm^3]
			#
			# where sigmaSB = 5.670373E-5 erg/cm^2/s/K^4. Then,
			#
			# Omega_gamma = rho_gamma / (Msun/g) * (kpc/cm)^3 / h^2 / AST_rho_crit_0_kpc3
			#
			# Most of these steps can be summarized in one constant.
			self.Ogamma0 = 4.48131796342E-07 * self.Tcmb0**4 / self.h2
			
			# The energy density in neutrinos is 7/8 (4/11)^(4/3) times the energy density in 
			# photons, per effective neutrino species.
			self.Onu0 = 0.22710731766 * self.Neff * self.Ogamma0
			
			# The density of relativistic species is the sum of the photon and neutrino densities.
			self.Or0 = self.Ogamma0 + self.Onu0
			
			# For convenience, compute the epoch of matter-radiation equality
			self.a_eq = self.Or0 / self.Om0
		else:
			self.Ogamma0 = 0.0
			self.Onu0 = 0.0
			self.Or0 = 0.0

		# Make sure flatness is obeyed
		self._ensureConsistency()
		
		# Flag for interpolation tables, storage, printing etc
		self.interpolation = interpolation
		self.storage = storage
		self.text_output = text_output
		self.print_info = print_info
		self.print_warnings = print_warnings
		
		# Lookup table for functions of z. This table runs from the future (a = 200.0) to 
		# a = 0.005. Due to some interpolation errors at the extrema of the range, the table 
		# runs to slightly lower and higher z than the interpolation is allowed for.
		self.z_min = -0.995
		self.z_min_compute = -0.998
		self.z_max = 200.01
		self.z_max_compute = 500.0
		self.z_Nbins = 50
		
		# Lookup table for P(k). The Pk_norm field is only needed if interpolation == False.
		# Note that the binning is highly irregular for P(k), since much more resolution is
		# needed at the BAO scale and around the bend in the power spectrum. Thus, the binning
		# is split into multiple regions with different resolutions.
		self.k_Pk = [1E-20, 1E-4, 5E-2, 1E0, 1E6, 1E20]
		self.k_Pk_Nbins = [10, 30, 60, 20, 10]
		
		# Lookup table for sigma. Note that the nominal accuracy to which the integral is 
		# evaluated should match with the accuracy of the interpolation which is set by Nbins.
		# Here, they are matched to be accurate to better than ~3E-3.
		self.R_min_sigma = 1E-12
		self.R_max_sigma = 1E3
		self.R_Nbins_sigma = 18.0
		self.accuracy_sigma = 3E-3
	
		# Lookup table for correlation function xi
		self.R_xi = [1E-3, 5E1, 5E2]
		self.R_xi_Nbins = [30, 40]
		self.accuracy_xi = 1E-5
		
		# Some functions permanently store lookup tables for faster execution. For this purpose,
		# we compute and save the storage path. 
		if self.storage:
			self.cache_dir = Utilities.getCacheDir(module = 'cosmology')
		
		# Note that the storage is active even if interpolation == False or storage == False, 
		# since fields can still be stored non-persistently (without writing to file).
		self._resetStorage()

		return

	###############################################################################################

	def __str__(self):
		
		s = 'Cosmology "%s", flat = %s, relspecies = %s, \n' \
			'    Om0 = %.4f, OL0 = %.4f, Ob0 = %.4f, H0 = %.2f, sigma8 = %.4f, ns = %.4f, \n' \
			'    Tcmb0 = %.4f, Neff = %.4f, PL = %s, PLn = %.4f' \
			% (self.name, str(self.flat), str(self.relspecies), \
			self.Om0, self.OL0, self.Ob0, self.H0, self.sigma8, self.ns, self.Tcmb0, self.Neff, \
			str(self.power_law), self.power_law_n)
		
		return s

	###############################################################################################

	def checkForChangedCosmology(self):
		"""
		Check whether the cosmological parameters have been changed by the user. If there are 
		changes, all pre-computed quantities (e.g., interpolation tables) are discarded and 
		re-computed if necessary.
		"""
		
		hash_new = self._getHash()
		if hash_new != self.hash_current:
			if self.print_warnings:
				print("Cosmology: Detected change in cosmological parameters.")
			self._ensureConsistency()
			self._resetStorage()
			
		return
	
	###############################################################################################
	# Utilities for internal use
	###############################################################################################
	
	# Depending on whether the cosmology is flat or not, OL0 and Ok0 take on certain values.

	def _ensureConsistency(self):
		
		if self.flat:
			self.OL0 = 1.0 - self.Om0 - self.Or0
			self.Ok0 = 0.0
		else:
			self.Ok0 = 1.0 - self.OL0 - self.Om0 - self.Or0

		return

	###############################################################################################

	# Compute a unique hash for the current cosmology name and parameters. If any of them change,
	# the hash will change, causing an update of stored quantities.
		
	def _getHash(self):
	
		param_string = "Name_%s_Flat_%s_relspecies_%s_Om0_%.4f_OL0_%.4f_Ob0_%.4f_H0_%.4f_sigma8_%.4f_ns_%.4f_Tcmb0_%.4f_Neff_%.4f_PL_%s_PLn_%.4f" \
			% (self.name, str(self.flat), str(self.relspecies), \
			self.Om0, self.OL0, self.Ob0, self.H0, self.sigma8, self.ns, self.Tcmb0, self.Neff, \
			str(self.power_law), self.power_law_n)

		hash_new = hashlib.md5(param_string).hexdigest()
	
		return hash_new
	
	###############################################################################################

	# Create a file name that is unique to this cosmology. While the hash encodes all necessary
	# information, the cosmology name is added to make it easier to identify the files with a 
	# cosmology.

	def _getUniqueFilename(self):
		
		return self.cache_dir + self.name + '_' + self._getHash()
	
	###############################################################################################

	# Load stored objects. This function is called during the __init__() routine, and if a change
	# in cosmological parameters is detected.

	def _resetStorage(self):

		# Reset the test hash and storage containers. There are two containes, one for objects
		# that are stored in a pickle file, and one for those that will be discarded when the 
		# class is destroyed.
		self.hash_current = self._getHash()
		self.storage_pers = {}
		self.storage_temp = {}
		
		# Check if there is a persistent object storage file. If so, load its contents into the
		# storage dictionary. We only load from file if the user has not switched of storage, and
		# if the user has not switched off interpolation.
		if self.storage and self.interpolation:
			filename_pickle = self._getUniqueFilename()
			if os.path.exists(filename_pickle):
				input_file = open(filename_pickle, "rb")
				self.storage_pers = pickle.load(input_file)
				input_file.close()

		return
	
	###############################################################################################

	# Permanent storage system for objects such as 2-dimensional data tables. If an object is 
	# already stored in memory, return it. If not, try to load it from file, otherwise return None.
	# Certain operations can already be performed on certain objects, so that they do not need to 
	# be repeated unnecessarily, for example:
	#
	# interpolator = True	Instead of a 2-dimensional table, return a spline interpolator that can
	#                       be used to evaluate the table.
	# inverse = True        Return an interpolator that gives x(y) instead of y(x)
	
	def _getStoredObject(self, object_name, interpolator = False, inverse = False):
		
		# Check for cosmology change
		self.checkForChangedCosmology()

		# Compute object name
		object_id = object_name
		if interpolator:
			object_id += '_interpolator'
		if inverse:
			object_id += '_inverse'

		# Find the object. There are multiple possibilities:
		# - Check for the exact object the user requested (the object_id)
		#   - Check in persistent storage
		#   - Check in temporary storage (where interpolator / inverse objects live)
		#   - Check in user text files
		# - Check for the raw object (the object_name)
		#   - Check in persistent storage
		#   - Check in user text files
		#   - Convert to the exact object, store in temporary storage
		# - If all fail, return None

		if object_id in self.storage_pers:	
			object_data = self.storage_pers[object_id]
		
		elif object_id in self.storage_temp:	
			object_data = self.storage_temp[object_id]

		elif self.storage and os.path.exists(self.cache_dir + object_id):
			object_data = numpy.loadtxt(self.cache_dir + object_id, usecols = (0, 1), \
									skiprows = 0, unpack = True)
			self.storage_temp[object_id] = object_data
			
		else:

			# We could not find the object ID anywhere. This can have two reasons: the object does
			# not exist, or we must transform an existing object.
			
			if interpolator:
				
				# First, a safety check; no interpolation objects should ever be requested if 
				# the user has switched off interpolation.
				if not self.interpolation:
					raise Exception('An interpolator object was requested even though interpolation is off.')
				
				# Try to find the object to transform. This object CANNOT be in temporary storage,
				# but it can be in persistent or user storage.
				object_raw = None
				
				if object_name in self.storage_pers:	
					object_raw = self.storage_pers[object_name]
		
				elif self.storage and os.path.exists(self.cache_dir + object_name):
					object_raw = numpy.loadtxt(self.cache_dir + object_name, usecols = (0, 1), \
									skiprows = 0, unpack = True)

				if object_raw is None:
					
					# We cannot find an object to convert, return none.
					object_data = None
				
				else:
					
					# Convert and store in temporary storage.
					if inverse: 
						
						# There is a subtlety: the spline interpolator can't deal with decreasing 
						# x-values, so if the y-values are decreasing, we reverse their order.
						if object_raw[1][-1] < object_raw[1][0]:
							object_raw = object_raw[:,::-1]
						
						object_data = scipy.interpolate.InterpolatedUnivariateSpline(object_raw[1], \
																					object_raw[0])
					else:
						object_data = scipy.interpolate.InterpolatedUnivariateSpline(object_raw[0], \
																					object_raw[1])
					self.storage_temp[object_id] = object_data
						
			else:
							
				# The object is not in storage at all, and cannot be generated; return none.
				object_data = None
				
		return object_data
	
	###############################################################################################

	# Save an object in memory and file storage. If persistent == True, this object is written to 
	# file storage (unless storage == False), and will be loaded the next time the same cosmology
	# is loaded. If persistent == False, the object is stored non-persistently.
	#
	# Note that all objects are reset if the cosmology changes. Thus, this function should be used
	# for ALL data that depend on cosmological parameters.
	
	def _storeObject(self, object_name, object_data, persistent = True):

		if persistent:
			self.storage_pers[object_name] = object_data
			
			if self.storage:
				# If the user has chosen text output, write a text file.
				if self.text_output:
					filename_text =  self.cache_dir + object_name
					numpy.savetxt(filename_text, numpy.transpose(object_data), fmt = "%.8e")
			
				# Store in file. We do not wish to save the entire storage dictionary, as there might be
				# user-defined objects in it.
				filename_pickle = self._getUniqueFilename()
				output_file = open(filename_pickle, "wb")
				pickle.dump(self.storage_pers, output_file, pickle.HIGHEST_PROTOCOL)
				output_file.close()  

		else:
			self.storage_temp[object_name] = object_data

		return
	
	###############################################################################################
	# Basic cosmology calculations
	###############################################################################################
	
	def Ez(self, z):
		"""
		The Hubble parameter as a function of redshift, in units of :math:`H_0`.
		
		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		E: array_like
			:math:`H(z) / H_0`; has the same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		Hz: The Hubble parameter as a function of redshift.
		"""
		
		ai = (1.0 + z)
		sum = self.Om0 * ai**3 + self.OL0
		if not self.flat:
			sum += self.Ok0 * ai**2
		if self.relspecies:
			sum += self.Or0 * ai**4
		E = numpy.sqrt(sum)
		
		return E

	###############################################################################################
	
	def Hz(self, z):
		"""
		The Hubble parameter as a function of redshift.
		
		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		H: array_like
			:math:`H(z)` in units of km/s/Mpc; has the same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		Ez: The Hubble parameter as a function of redshift, in units of :math:`H_0`.
		"""
		
		H = self.Ez(z) * self.H0
					
		return H

	###############################################################################################

	# Standard cosmological integrals. These integrals are not persistently stored in files because
	# they can be evaluated between any two redshifts which would make the tables very large.
	#
	# z_min and z_max can be numpy arrays or numbers. If one of the two is a number and the other an
	# array, the same z_min / z_max is used for all z_min / z_max in the array (this is useful if 
	# z_max = inf, for example).
	
	def _integral(self, integrand, z_min, z_max):

		min_is_array = Utilities.isArray(z_min)
		max_is_array = Utilities.isArray(z_max)
		use_array = min_is_array or max_is_array
		
		if use_array and not min_is_array:
			z_min_use = numpy.array([z_min] * len(z_max))
		else:
			z_min_use = z_min
		
		if use_array and not max_is_array:
			z_max_use = numpy.array([z_max] * len(z_min))
		else:
			z_max_use = z_max
		
		if use_array:
			if min_is_array and max_is_array and len(z_min) != len(z_max):
				raise Exception("If both z_min and z_max are arrays, they need to have the same size.")
			integ = numpy.zeros((len(z_min_use)), numpy.float)
			for i in range(len(z_min_use)):
				integ[i], _ = scipy.integrate.quad(integrand, z_min_use[i], z_max_use[i])
		else:
			integ, _ = scipy.integrate.quad(integrand, z_min, z_max)
		
		return integ
	
	###############################################################################################

	# The integral over 1 / E(z) enters into the comoving distance.

	def _integral_oneOverEz(self, z_min, z_max = numpy.inf):
		
		def integrand(z):
			return 1.0 / self.Ez(z)
		
		return self._integral(integrand, z_min, z_max)

	###############################################################################################

	# The integral over 1 / E(z) / (1 + z) enters into the age of the universe.

	def _integral_oneOverEz1pz(self, z_min, z_max = numpy.inf):
		
		def integrand(z):
			return 1.0 / self.Ez(z) / (1.0 + z)
		
		return self._integral(integrand, z_min, z_max)

	###############################################################################################

	# Used by _zFunction

	def _zInterpolator(self, table_name, func, inverse = False, future = True):

		table_name = table_name + '_%s' % (self.name) 
		interpolator = self._getStoredObject(table_name, interpolator = True, inverse = inverse)
		
		if interpolator is None:
			if self.print_info:
				print("Computing lookup table in z.")
			
			if future:
				log_min = numpy.log10(1.0 + self.z_min_compute)
			else:
				log_min = 0.0
			log_max = numpy.log10(1.0 + self.z_max_compute)
			bin_width = (log_max - log_min) / self.z_Nbins
			z_table = 10**numpy.arange(log_min, log_max + bin_width, bin_width) - 1.0
			x_table = func(z_table)
			
			self._storeObject(table_name, numpy.array([z_table, x_table]))
			if self.print_info:
				print("Lookup table completed.")
			interpolator = self._getStoredObject(table_name, interpolator = True, inverse = inverse)
		
		return interpolator

	###############################################################################################

	# General container for methods that are functions of z and use interpolation
	
	def _zFunction(self, table_name, func, z, inverse = False, future = True, derivative = 0):

		if self.interpolation:
			
			# Get interpolator. If it does not exist, create it.
			interpolator = self._zInterpolator(table_name, func, inverse = inverse, future = future)
			
			# Check limits of z array. If inverse == True, we need to check the limits on 
			# the result function.
			if inverse:
				min_ = interpolator.get_knots()[0]
				max_ = interpolator.get_knots()[-1]
				
				if numpy.min(z) < min_:
					msg = "f = %.3f outside range (min. f is %.3f)." \
						% (numpy.min(z), min_)
					raise Exception(msg)
					
				if numpy.max(z) > max_:
					msg = "f = %.3f outside range (max. f is %.3f)." \
						% (numpy.max(z), max_)
					raise Exception(msg)
				
			else:
				if numpy.min(z) < self.z_min:
					msg = "z = %.2f outside range (min. z is %.2f)." \
						% (numpy.min(z), self.z_min)
					raise Exception(msg)
					
				if numpy.max(z) > self.z_max:
					msg = "z = %.2f outside range (max. z is %.2f)." \
						% (numpy.max(z), self.z_max)
					raise Exception(msg)
				
			ret = interpolator(z, nu = derivative)				
			
		else:
			if derivative > 0:
				raise Exception("Derivative can only be evaluated if interpolation == True.")

			ret = func(z)
		
		return ret

	###############################################################################################
	# Times & distances
	###############################################################################################
	
	def hubbleTime(self, z = 0.0):
		"""
		The Hubble time, :math:`1/H(z)`.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		tH: float
			:math:`1/H` in units of Gyr; has the same dimensions as z. By default returns 
			:math:`1/H_0`.

		See also
		-------------------------------------------------------------------------------------------
		lookbackTime: The lookback time since z.
		age: The age of the universe at redshift z.
		"""
		
		tH = 1E-16 * AST_Mpc / AST_year / self.h / self.Ez(z)
		
		return tH
	
	###############################################################################################

	def _lookbackTimeExact(self, z):
		
		t = self.hubbleTime() * self._integral_oneOverEz1pz(0.0, z)

		return t

	###############################################################################################

	def lookbackTime(self, z, derivative = 0, inverse = False):
		"""
		The lookback time since z.
		
		The lookback time corresponds to the difference between the age of the universe at 
		redshift z and today.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift, where :math:`-0.995 < z < 200`; can be a number or a numpy array.
		derivative: int
			If greater than 0, evaluate the nth derivative, :math:`d^nt/dz^n`.
		inverse: bool
			If True, evaluate :math:`z(t)` instead of :math:`t(z)`.

		Returns
		-------------------------------------------------------------------------------------------
		t: array_like
			The lookback time (or its derivative) since z in units of Gigayears; has the same 
			dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		hubbleTime: The Hubble time, :math:`1/H_0`.
		age: The age of the universe at redshift z.
		"""
		
		t = self._zFunction('lookbacktime', self._lookbackTimeExact, z, derivative = derivative, \
						inverse = inverse)
		
		return t
	
	###############################################################################################

	def _ageExact(self, z):
		
		t = self.hubbleTime() * self._integral_oneOverEz1pz(z, numpy.inf)
		
		return t
	
	###############################################################################################
	
	def age(self, z = 0.0, derivative = 0, inverse = False):
		"""
		The age of the universe at redshift z.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift, where :math:`-0.995 < z < 200`; can be a number or a numpy array.
		derivative: int
			If greater than 0, evaluate the nth derivative, :math:`d^nt/dz^n`.
		inverse: bool
			If True, evaluate :math:`z(t)` instead of :math:`t(z)`.

		Returns
		-------------------------------------------------------------------------------------------
		t: array_like
			The age of the universe (or its derivative) at redshift z in Gigayears; has the 
			same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		hubbleTime: The Hubble time, :math:`1/H_0`.
		lookbackTime: The lookback time since z.
		"""

		t = self._zFunction('age', self._ageExact, z, derivative = derivative, inverse = inverse)
		
		return t
	
	###############################################################################################
	
	# This function does not use interpolation because both zmin and zmax are free, which would 
	# lead to a more complicated 2D-interpolation.
	
	def comovingDistance(self, z_min = 0.0, z_max = 0.0):
		"""
		The comoving distance between redshift :math:`z_{min}` and :math:`z_{max}`.
		
		Either z_min or z_min can be a numpy array; in those cases, the same z_min / z_max is 
		applied to all values of the other. If both are numpy arrays, they need to have 
		the same dimensions, and the comoving distance returned corresponds to a series of 
		different z_min and z_max values. 

		Parameters
		-------------------------------------------------------------------------------------------
		zmin: array_like
			Redshift; can be a number or a numpy array.
		zmax: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		d: array_like
			The comoving distance in Mpc/h; has the same dimensions as zmin and/or zmax.

		See also
		-------------------------------------------------------------------------------------------
		luminosityDistance: The luminosity distance to redshift z.
		angularDiameterDistance: The angular diameter distance to redshift z.
		"""
		
		d = self._integral_oneOverEz(z_min = z_min, z_max = z_max) * AST_c * 1E-7
		
		return d

	###############################################################################################

	def _luminosityDistanceExact(self, z):
		
		d = self.comovingDistance(z_min = 0.0, z_max = z) * (1.0 + z)
		
		return d

	###############################################################################################
	
	def luminosityDistance(self, z, derivative = 0, inverse = False):
		"""
		The luminosity distance to redshift z.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift, where :math:`-0.995 < z < 200`; can be a number or a numpy array.
		derivative: int
			If greater than 0, evaluate the nth derivative, :math:`d^nD/dz^n`.
		inverse: bool
			If True, evaluate :math:`z(D)` instead of :math:`D(z)`.
			
		Returns
		-------------------------------------------------------------------------------------------
		d: array_like
			The luminosity distance (or its derivative) in Mpc/h; has the same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		comovingDistance: The comoving distance between redshift :math:`z_{min}` and :math:`z_{max}`.
		angularDiameterDistance: The angular diameter distance to redshift z.
		"""
		
		d = self._zFunction('luminositydist', self._luminosityDistanceExact, z, \
						future = False, derivative = derivative, inverse = inverse)
		
		return d
	
	###############################################################################################

	def _angularDiameterDistanceExact(self, z):
		
		d = self.comovingDistance(z_min = 0.0, z_max = z) / (1.0 + z)
		
		return d

	###############################################################################################

	def angularDiameterDistance(self, z, derivative = 0, inverse = False):
		"""
		The angular diameter distance to redshift z.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift, where :math:`-0.995 < z < 200`; can be a number or a numpy array.
		derivative: int
			If greater than 0, evaluate the nth derivative, :math:`d^nD/dz^n`.
		inverse: bool
			If True, evaluate :math:`z(D)` instead of :math:`D(z)`.

		Returns
		-------------------------------------------------------------------------------------------
		d: array_like
			The angular diameter distance (or its derivative) in Mpc/h; has the same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		comovingDistance: The comoving distance between redshift :math:`z_{min}` and :math:`z_{max}`.
		luminosityDistance: The luminosity distance to redshift z.
		"""

		d = self._zFunction('angdiamdist', self._angularDiameterDistanceExact, z, \
						future = False, derivative = derivative, inverse = inverse)
		
		return d

	###############################################################################################

	# This function is not interpolated because the distance modulus is not defined at z = 0.

	def distanceModulus(self, z):
		"""
		The distance modulus to redshift z in magnitudes.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		mu: array_like
			The distance modulus in magnitudes; has the same dimensions as z.
		"""
		
		mu = 5.0 * numpy.log10(self.luminosityDistance(z) / self.h * 1E5)
		
		return mu

	###############################################################################################

	def soundHorizon(self):
		"""
		The sound horizon at recombination.

		This function returns the sound horizon in Mpc (not Mpc/h!), according to Eisenstein & Hu 
		1998, equation 26. This fitting function is accurate to 2% where :math:`\Omega_b h^2 > 0.0125` 
		and :math:`0.025 < \Omega_m h^2 < 0.5`.

		Returns
		-------------------------------------------------------------------------------------------
		s: float
			The sound horizon at recombination in Mpc.
		"""
				
		s = 44.5 * numpy.log(9.83 / self.Omh2) / numpy.sqrt(1.0 + 10.0 * self.Ombh2**0.75)
		
		return s

	###############################################################################################
	# Densities and overdensities
	###############################################################################################
	
	def rho_c(self, z):
		"""
		The critical density of the universe at redshift z.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		rho_critical: array_like
			The critical density in units of physical :math:`M_{\odot} h^2 / kpc^3`; has the same 
			dimensions as z.
		"""
			
		return AST_rho_crit_0_kpc3 * self.Ez(z)**2

	###############################################################################################
	
	def rho_m(self, z):
		"""
		The matter density of the universe at redshift z.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		rho_matter: array_like
			The matter density in units of physical :math:`M_{\odot} h^2 / kpc^3`; has the same 
			dimensions as z.
	
		See also
		-------------------------------------------------------------------------------------------
		Om: The matter density of the universe, in units of the critical density.
		"""
			
		return AST_rho_crit_0_kpc3 * self.Om0 * (1.0 + z)**3


	###############################################################################################
	
	def rho_L(self, z):
		"""
		The dark energy density of the universe at redshift z.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		rho_Lambda: array_like
			The dark energy density in units of physical :math:`M_{\odot} h^2 / kpc^3`; has the same 
			dimensions as z.
	
		See also
		-------------------------------------------------------------------------------------------
		OL: The dark energy density of the universe, in units of the critical density. 
		"""
			
		return AST_rho_crit_0_kpc3 * self.OL0

	###############################################################################################
	
	def rho_gamma(self, z):
		"""
		The photon density of the universe at redshift z.
		
		If ``relspecies == False``, this function returns 0.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		rho_gamma: array_like
			The photon density in units of physical :math:`M_{\odot} h^2 / kpc^3`; has the same 
			dimensions as z.
	
		See also
		-------------------------------------------------------------------------------------------
		Ogamma: The density of photons in the universe, in units of the critical density.
		"""
			
		return AST_rho_crit_0_kpc3 * self.Ogamma0 * (1.0 + z)**4

	###############################################################################################
	
	def rho_nu(self, z):
		"""
		The neutrino density of the universe at redshift z.

		If ``relspecies == False``, this function returns 0.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		rho_nu: array_like
			The neutrino density in units of physical :math:`M_{\odot} h^2 / kpc^3`; has the same 
			dimensions as z.
	
		See also
		-------------------------------------------------------------------------------------------
		Onu: The density of neutrinos in the universe, in units of the critical density.
		"""
			
		return AST_rho_crit_0_kpc3 * self.Onu0 * (1.0 + z)**4

	###############################################################################################
	
	def rho_r(self, z):
		"""
		The density of relativistic species in the universe at redshift z.
		
		This density is the sum of the photon and neutrino densities. If ``relspecies == False``, 
		this function returns 0.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		rho_relativistic: array_like
			The density of relativistic species in units of physical :math:`M_{\odot} h^2 / kpc^3`; 
			has the same dimensions as z.
	
		See also
		-------------------------------------------------------------------------------------------
		Or: The density of relativistic species in the universe, in units of the critical density.
		"""
			
		return AST_rho_crit_0_kpc3 * self.Or0 * (1.0 + z)**4

	###############################################################################################

	def Om(self, z):
		"""
		The matter density of the universe, in units of the critical density.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		Omega_matter: array_like
			Has the same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		rho_m: The matter density of the universe at redshift z.
		"""

		return self.Om0 * (1.0 + z)**3 / (self.Ez(z))**2

	###############################################################################################

	def OL(self, z):
		"""
		The dark energy density of the universe, in units of the critical density. 
		
		In a flat universe, :math:`\Omega_{\Lambda} = 1 - \Omega_m - \Omega_r`.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		Omega_Lambda: array_like
			Has the same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		rho_L: The dark energy density of the universe at redshift z.
		"""

		return self.OL0 / (self.Ez(z))**2

	###############################################################################################

	def Ok(self, z):
		"""
		The curvature density of the universe in units of the critical density. 
		
		In a flat universe, :math:`\Omega_k = 0`.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		Omega_curvature: array_like
			Has the same dimensions as z.
		"""
					
		return self.Ok0 * (1.0 + z)**2 / (self.Ez(z))**2

	###############################################################################################

	def Ogamma(self, z):
		"""
		The density of photons in the universe, in units of the critical density.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		Omega_gamma: array_like
			Has the same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		rho_gamma: The photon density of the universe at redshift z.
		"""
					
		return self.Ogamma0 * (1.0 + z)**4 / (self.Ez(z))**2

	###############################################################################################

	def Onu(self, z):
		"""
		The density of neutrinos in the universe, in units of the critical density.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		Omega_nu: array_like
			Has the same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		rho_nu: The neutrino density of the universe at redshift z.
		"""
					
		return self.Onu0 * (1.0 + z)**4 / (self.Ez(z))**2
	
	###############################################################################################

	def Or(self, z):
		"""
		The density of relativistic species in the universe, in units of the critical density. 

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		Omega_relativistic: array_like
			Has the same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		rho_r: The density of relativistic species in the universe at redshift z.
		"""
					
		return self.Or0 * (1.0 + z)**4 / (self.Ez(z))**2

	###############################################################################################
	# Structure growth, power spectrum etc.
	###############################################################################################
	
	def lagrangianR(self, M):
		"""
		The lagrangian radius of a halo of mass M.

		Converts the mass of a halo (in comoving :math:`M_{\odot} / h`) to the radius of its 
		comoving Lagrangian volume (in comoving Mpc/h), that is the volume that encloses the halo's 
		mass given the mean density of the universe at z = 0.

		Parameters
		-------------------------------------------------------------------------------------------
		M: array_like
			Halo mass in :math:`M_{\odot} / h`; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		R: array_like
			The lagrangian radius in comoving Mpc/h; has the same dimensions as M.

		See also
		-------------------------------------------------------------------------------------------
		lagrangianM: The lagrangian mass of a halo of radius R.
		"""
		
		return (3.0 * M / 4.0 / math.pi / self.rho_m(0.0) / 1E9)**(1.0 / 3.0)
	
	###############################################################################################
	
	def lagrangianM(self, R):
		"""
		The lagrangian mass of a halo of radius R.

		Converts the radius of a halo (in comoving Mpc/h) to the mass in its comoving Lagrangian 
		volume (in :math:`M_{\odot} / h`), that is the volume that encloses the halo's mass given the 
		mean density of the universe at z = 0.

		Parameters
		-------------------------------------------------------------------------------------------
		R: array_like
			Halo radius in comoving Mpc/h; can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		M: array_like
			The lagrangian mass; has the same dimensions as R.

		See also
		-------------------------------------------------------------------------------------------
		lagrangianR: The lagrangian radius of a halo of mass M.
		"""
				
		return 4.0 / 3.0 * math.pi * R**3 * self.rho_m(0.0) * 1E9

	###############################################################################################

	def growthFactorUnnormalized(self, z):
		"""
		The linear growth factor, :math:`D_+(z)`.
		
		The growth factor describes the linear evolution of over- and underdensities in the dark
		matter density field. There are three regimes: 1) In the matter-radiation regime, we use an 
		approximate analytical formula (Equation 5 in Gnedin, Kravtsov & Rudd 2011). If relativistic 
		species are ignored, :math:`D_+(z) \propto a`. 2) In the matter-dominated regime, 
		:math:`D_+(z) \propto a`. 3) In the matter-dark energy regime, we evaluate :math:`D_+(z)` 
		through integration as defined in Eisenstein & Hu 99, Equation 8 (see also Heath 1977). 
		
		At the transition between the integral and analytic approximation regimes, the two 
		expressions do not quite match up, with differences of the order <1E-3. in order to avoid
		a discontinuity, we introduce a transition regime where the two quantities are linearly
		interpolated.
		
		The normalization is such that the growth factor approaches :math:`D_+(a) = a` in the 
		matter-dominated regime. There are other normalizations of the growth factor (e.g., Percival 
		2005, Equation 15), but since we almost always care about the growth factor normalized to 
		z = 0, the normalization does not matter too much (see the :func:`growthFactor` function).
		
		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift, where :math:`-0.995 < z`; the high end of z is only limited by the validity 
			of the analytical approximation mentioned above. Can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		D: array_like
			The linear growth factor; has the same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		growthFactor: The linear growth factor normalized to z = 0, :math:`D_+(z) / D_+(0)`.

		Warnings
		-------------------------------------------------------------------------------------------
		This function directly evaluates the growth factor by integration or analytical 
		approximation. In most cases, the :func:`growthFactor` function should be used since it 
		interpolates and is thus much faster.
		"""

		# The growth factor integral uses E(z), but is not designed to take relativistic species
		# into account. Thus, using the standard E(z) leads to wrong results. Instead, we pretend
		# that the small radiation content at low z behaves like dark energy which leads to a very
		# small error but means that the formula converges to a at high z.
		def Ez_D(z):
			ai = (1.0 + z)
			sum = self.Om0 * ai**3 + self.OL0
			if self.relspecies:
				sum += self.Or0
			if not self.flat:
				sum += self.Ok0 * ai**2
			E = numpy.sqrt(sum)
			return E

		# The integrand
		def integrand(z):
			return (1.0 + z) / (Ez_D(z))**3

		# Create a transition regime centered around z = 10 in log space
		z_switch = 10.0
		trans_width = 2.0
		zt1 = z_switch * trans_width
		zt2 = z_switch / trans_width
		
		# Split into late (1), early (2) and a transition interval (3)
		z_arr, is_array = Utilities.getArray(z)
		a = 1.0 / (1.0 + z_arr)
		D = numpy.zeros((len(z_arr)), numpy.float)
		mask1 = z_arr < (zt1)
		mask2 = z_arr > (zt2)
		mask3 = mask1 & mask2
		
		# Compute D from integration at low redshift
		z1 = z_arr[mask1]
		D[mask1] = 5.0 / 2.0 * self.Om0 * Ez_D(z1) * self._integral(integrand, z1, numpy.inf)
		D1 = D[mask3]
		
		# Compute D analytically at high redshift.
		a2 = a[mask2]
		if self.relspecies:
			x = a2 / self.a_eq
			term1 = numpy.sqrt(1.0 + x)
			term2 = 2.0 * term1 + (2.0 / 3.0 + x) * numpy.log((term1 - 1.0) / (term1 + 1.0))
			D[mask2] = a2 + 2.0 / 3.0 * self.a_eq + self.a_eq / (2.0 * numpy.log(2.0) - 3.0) * term2
		else:
			D[mask2] = a2
		D2 = D[mask3]

		# Average in transition regime
		at1 = numpy.log(1.0 / (zt1 + 1.0))
		at2 = numpy.log(1.0 / (zt2 + 1.0))
		dloga = at2 - at1
		loga = numpy.log(a[mask3])
		D[mask3] = (D1 * (loga - at1) + D2 * (at2 - loga)) / dloga

		# Reduce array to number if necessary
		if not is_array:
			D = D[0]
		
		return D

	###############################################################################################

	def _growthFactorExact(self, z):
		
		D = self.growthFactorUnnormalized(z) / self.growthFactorUnnormalized(0.0)
		
		return D

	###############################################################################################

	def growthFactor(self, z, derivative = 0, inverse = False):
		"""
		The linear growth factor normalized to z = 0, :math:`D_+(z) / D_+(0)`.

		The growth factor describes the linear evolution of over- and underdensities in the dark
		matter density field. This function is sped up through interpolation which barely degrades 
		its accuracy, but if you wish to evaluate the exact integral or compute the growth factor 
		for very high redshifts (z > 200), please use the :func:`growthFactorUnnormalized`
		function.

		Parameters
		-------------------------------------------------------------------------------------------
		z: array_like
			Redshift, where :math:`-0.995 < z < 200`; can be a number or a numpy array.
		derivative: int
			If greater than 0, evaluate the nth derivative, :math:`d^nD_+/dz^n`.
		inverse: bool
			If True, evaluate :math:`z(D_+)` instead of :math:`D_+(z)`.

		Returns
		-------------------------------------------------------------------------------------------
		D: array_like
			The linear growth factor (or its derivative); has the same dimensions as z.

		See also
		-------------------------------------------------------------------------------------------
		growthFactorUnnormalized: The linear growth factor, :math:`D_+(z)`.
		"""

		D = self._zFunction('growthfactor', self._growthFactorExact, z, derivative = derivative, \
						inverse = inverse)

		return D

	###############################################################################################
	
	def collapseOverdensity(self, deltac_const = True, sigma = None):
		"""
		The threshold overdensity for halo collapse.
		
		For most applications, ``deltac_const = True`` works fine; in that case, this function
		simply returns the collapse overdensity predicted by the top-hat collapse model, 1.686. 
		Alternatively, a correction for the ellipticity of peaks can be applied according to Sheth 
		et al. 2001. In that case, the variance on the scale of a halo must also be passed.

		Parameters
		-------------------------------------------------------------------------------------------
		deltac_const: bool
			If True, the function returns the constant top-hat model collapse overdensity. If False,
			a correction due to the ellipticity of halos is applied.
		sigma: float
			The rms variance on the scale of the halo; only necessary if ``deltac_const == False``.

		Returns
		-------------------------------------------------------------------------------------------
		delta_c: float
			The threshold overdensity for collapse.

		See also
		-------------------------------------------------------------------------------------------
		sigma: The rms variance of the linear density field on a scale R, :math:`\\sigma(R)`.
		"""
				
		if deltac_const:
			delta_c = AST_delta_collapse
		else:
			delta_c = AST_delta_collapse * (1.0 + 0.47 * (sigma / AST_delta_collapse)**1.23)
		
		return delta_c
	
	###############################################################################################
	
	def transferFunctionEH98(self, k):
		"""
		The transfer function according to Eisenstein & Hu 1998.
		
		The transfer function transforms the spectrum of primordial fluctuations into the
		power spectrum of the initial matter density fluctuations. The primordial power spectrum is 
		usually described as a power law, leading to a power spectrum
		
		.. math::
			P(k) = T(k)^2 k^{n_s}
			
		where P(k) is the matter power spectrum, T(k) is the transfer function, and :math:`n_s` is 
		the tilt of the primordial power spectrum. This function computes the Eisenstein & Hu 1998 
		approximation to the transfer function at a scale k, and is based on Matt Becker's 
		cosmocalc code.
	
		Parameters
		-------------------------------------------------------------------------------------------
		k: array_like
			The wavenumber k (in comoving h/Mpc); can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		Tk: array_like
			The transfer function; has the same dimensions as k.

		See also
		-------------------------------------------------------------------------------------------
		transferFunctionEH98Smooth: The transfer function according to Eisenstein & Hu 1998, without the BAO features.
		"""
					
		# Define shorter expressions
		omb = self.Ob0
		om0 = self.Om0
		omc = om0 - omb
		h = self.h
		theta2p7 = self.Tcmb0 / 2.7
		
		# Convert kh from h/Mpc to 1/Mpc
		kh = k * h
	
		# Equation 2
		zeq = 2.50e4 * om0 * h * h / (theta2p7 * theta2p7 * theta2p7 * theta2p7)
	
		# Equation 3
		keq = 7.46e-2 * om0 * h * h / (theta2p7 * theta2p7)
	
		# Equation 4
		b1d = 0.313 * pow(om0 * h * h, -0.419) * (1.0 + 0.607 * pow(om0 * h * h, 0.674))
		b2d = 0.238 * pow(om0 * h * h, 0.223)
		zd = 1291.0 * pow(om0 * h * h, 0.251) / (1.0 + 0.659 * pow(om0 * h * h, 0.828)) \
			* (1.0 + b1d * pow(omb * h * h, b2d))
	
		# Equation 5
		Rd = 31.5 * omb * h * h / (theta2p7 * theta2p7 * theta2p7 * theta2p7) / ((zd) / 1e3)
		Req = 31.5 * omb * h * h / (theta2p7 * theta2p7 * theta2p7 * theta2p7) / (zeq / 1e3)
	
		# Equation 6
		s = 2.0 / 3.0 / keq * numpy.sqrt(6.0 / Req) * numpy.log((numpy.sqrt(1.0 + Rd) + \
			numpy.sqrt(Rd + Req)) / (1.0 + numpy.sqrt(Req)))
	
		# Equation 7
		ksilk = 1.6 * pow(omb * h * h, 0.52) * pow(om0 * h * h, 0.73) \
			* (1.0 + pow(10.4 * om0 * h * h, -0.95))
	
		# Equation 10
		q = kh / 13.41 / keq
	
		# Equation 11
		a1 = pow(46.9 * om0 * h * h, 0.670) * (1.0 + pow(32.1 * om0 * h * h, -0.532))
		a2 = pow(12.0 * om0 * h * h, 0.424) * (1.0 + pow(45.0 * om0 * h * h, -0.582))
		ac = pow(a1, -1.0 * omb / om0) * pow(a2, -1.0 * (omb / om0) * (omb / om0) * (omb / om0))
	
		# Equation 12
		b1 = 0.944 / (1.0 + pow(458.0 * om0 * h * h, -0.708))
		b2 = pow(0.395 * om0 * h * h, -0.0266)
		bc = 1.0 / (1.0 + b1 * (pow(omc / om0, b2) - 1.0))
	
		# Equation 15
		y = (1.0 + zeq) / (1.0 + zd)
		Gy = y * (-6.0 * numpy.sqrt(1.0 + y) + (2.0 + 3.0 * y) \
			* numpy.log((numpy.sqrt(1.0 + y) + 1.0) / (numpy.sqrt(1.0 + y) - 1.0)))
	
		# Equation 14
		ab = 2.07 * keq * s * pow(1.0 + Rd, -3.0 / 4.0) * Gy
	
		# Get CDM part of transfer function
	
		# Equation 18
		f = 1.0 / (1.0 + (kh * s / 5.4) * (kh * s / 5.4) * (kh * s / 5.4) * (kh * s / 5.4))
	
		# Equation 20
		C = 14.2 / ac + 386.0 / (1.0 + 69.9 * pow(q, 1.08))
	
		# Equation 19
		T0t = numpy.log(numpy.e + 1.8 * bc * q) / (numpy.log(numpy.e + 1.8 * bc * q) + C * q * q)
	
		# Equation 17
		C1bc = 14.2 + 386.0 / (1.0 + 69.9 * pow(q, 1.08))
		T0t1bc = numpy.log(numpy.e + 1.8 * bc * q) / (numpy.log(numpy.e + 1.8 * bc * q) + C1bc * q * q)
		Tc = f * T0t1bc + (1.0 - f) * T0t
	
		# Get baryon part of transfer function
	
		# Equation 24
		bb = 0.5 + omb / om0 + (3.0 - 2.0 * omb / om0) * numpy.sqrt((17.2 * om0 * h * h) \
			* (17.2 * om0 * h * h) + 1.0)
	
		# Equation 23
		bnode = 8.41 * pow(om0 * h * h, 0.435)
	
		# Equation 22
		st = s / pow(1.0 + (bnode / kh / s) * (bnode / kh / s) * (bnode / kh / s), 1.0 / 3.0)
	
		# Equation 21
		C11 = 14.2 + 386.0 / (1.0 + 69.9 * pow(q, 1.08))
		T0t11 = numpy.log(numpy.e + 1.8 * q) / (numpy.log(numpy.e + 1.8 * q) + C11 * q * q)
		Tb = (T0t11 / (1.0 + (kh * s / 5.2) * (kh * s / 5.2)) + ab / (1.0 + (bb / kh / s) * \
			(bb / kh / s) * (bb / kh / s)) * numpy.exp(-pow(kh / ksilk, 1.4))) \
			* numpy.sin(kh * st) / (kh * st)
	
		# Total transfer function
		Tk = omb / om0 * Tb + omc / om0 * Tc
	
		return Tk

	###############################################################################################
	
	# The Eisenstein & Hu 1998 transfer function at a scale k (h / Mpc) but without the BAO wiggles.
	
	def transferFunctionEH98Smooth(self, k):
		"""
		The transfer function according to Eisenstein & Hu 1998, without the BAO features.
		
		Same as the :func:`transferFunctionEH98` function, but in a simplified version
		without the baryon acoustic oscillation (BAO) features.
	
		Parameters
		-------------------------------------------------------------------------------------------
		k: array_like
			The wavenumber k (in comoving h/Mpc); can be a number or a numpy array.

		Returns
		-------------------------------------------------------------------------------------------
		Tk: array_like
			The transfer function; has the same dimensions as k.

		See also
		-------------------------------------------------------------------------------------------
		transferFunctionEH98: The transfer function according to Eisenstein & Hu 1998.
		"""
		omb = self.Ob0
		om0 = self.Om0
		h = self.h
		theta2p7 = self.Tcmb0 / 2.7

		# Convert kh from hMpc^-1 to Mpc^-1
		kh = k * h
	
		# Equation 26
		s = 44.5 * numpy.log(9.83 / om0 / h / h) / numpy.sqrt(1.0 + 10.0 * pow(omb * h * h, 0.75))
	
		# Equation 31
		alphaGamma = 1.0 - 0.328 * numpy.log(431.0 * om0 * h * h) * omb / om0 \
				+ 0.38 * numpy.log(22.3 * om0 * h * h) * (omb / om0) * (omb / om0)
	
		# Equation 30
		Gamma = om0 * h * (alphaGamma + (1.0 - alphaGamma) / (1.0 + pow(0.43 * kh * s, 4.0)))
	
		# Equation 28
		q = k * theta2p7 * theta2p7 / Gamma
	
		# Equation 29
		C0 = 14.2 + 731.0 / (1.0 + 62.5 * q)
		L0 = numpy.log(2.0 * numpy.exp(1.0) + 1.8 * q)
		Tk = L0 / (L0 + C0 * q * q)
	
		return Tk	

	###############################################################################################

	def _matterPowerSpectrumExact(self, k, Pk_source = 'eh98', ignore_norm = False):

		if self.power_law:
			
			Pk_source = 'powerlaw'
			Pk = k**self.power_law_n
			
		elif Pk_source == 'eh98':

			T = self.transferFunctionEH98(k)
			Pk = T * T * k**self.ns

		elif Pk_source == 'eh98smooth':

			T = self.transferFunctionEH98Smooth(k)
			Pk = T * T * k**self.ns

		else:
			
			table_name = 'matterpower_%s_%s' % (self.name, Pk_source)
			table = self._getStoredObject(table_name)

			if table is None:
				msg = "Could not load data table, %s." % (table_name)
				raise Exception(msg)
			if numpy.max(k) > numpy.max(table[0]):
				msg = "k (%.2e) is larger than max. k (%.2e)." % (numpy.max(k), numpy.max(table[0]))
				raise Exception(msg)
			if numpy.min(k) < numpy.min(table[0]):
				msg = "k (%.2e) is smaller than min. k (%.2e)." % (numpy.min(k), numpy.min(table[0]))
				raise Exception(msg)

			Pk = numpy.interp(k, table[0], table[1])
		
		# This is a little tricky. We need to store the normalization factor somewhere, even if 
		# interpolation = False; otherwise, we get into an infinite loop of computing sigma8, P(k), 
		# sigma8 etc.
		if not ignore_norm:
			norm_name = 'Pk_norm_%s_%s' % (self.name, Pk_source)
			norm = self._getStoredObject(norm_name)
			if norm is None:
				sigma_8Mpc = self._sigmaExact(8.0, filt = 'tophat', Pk_source = Pk_source, \
											exact_Pk = True, ignore_norm = True)
				norm = (self.sigma8 / sigma_8Mpc)**2
				self._storeObject(norm_name, norm, persistent = False)

			Pk *= norm

		return Pk
	
	###############################################################################################

	# Utility to get the min and max k for which a power spectrum is valid. Only for internal use.

	def _matterPowerSpectrumLimits(self, Pk_source):
		
		if self.power_law or Pk_source == 'eh98' or Pk_source == 'eh98smooth':
			k_min = self.k_Pk[0]
			k_max = self.k_Pk[-1]
		else:
			table_name = 'matterpower_%s_%s' % (self.name, Pk_source)
			table = self._getStoredObject(table_name)

			if table is None:
				msg = "Could not load data table, %s." % (table_name)
				raise Exception(msg)
	
			k_min = table[0][0]
			k_max = table[0][-1]
		
		return k_min, k_max
	
	###############################################################################################

	# Return a spline interpolator for the power spectrum. Generally, P(k) should be evaluated 
	# using the matterPowerSpectrum() function below, but for some performance-critical operations
	# it is faster to obtain the interpolator directly from this function. Note that the lookup 
	# table created here is complicated, with extra resolution around the BAO scale.

	def _matterPowerSpectrumInterpolator(self, Pk_source, inverse = False):
		
		table_name = 'Pk_%s_%s' % (self.name, Pk_source)
		interpolator = self._getStoredObject(table_name, interpolator = True, inverse = inverse)
	
		if interpolator is None:
			if self.print_info:
				print("Cosmology.matterPowerSpectrum: Computing lookup table.")				
			data_k = numpy.zeros((numpy.sum(self.k_Pk_Nbins) + 1), numpy.float)
			n_regions = len(self.k_Pk_Nbins)
			k_computed = 0
			for i in range(n_regions):
				log_min = numpy.log10(self.k_Pk[i])
				log_max = numpy.log10(self.k_Pk[i + 1])
				log_range = log_max - log_min
				bin_width = log_range / self.k_Pk_Nbins[i]
				if i == n_regions - 1:
					data_k[k_computed:k_computed + self.k_Pk_Nbins[i] + 1] = \
						10**numpy.arange(log_min, log_max + bin_width, bin_width)
				else:
					data_k[k_computed:k_computed + self.k_Pk_Nbins[i]] = \
						10**numpy.arange(log_min, log_max, bin_width)
				k_computed += self.k_Pk_Nbins[i]
			
			data_Pk = self._matterPowerSpectrumExact(data_k, Pk_source = Pk_source, ignore_norm = False)
			table_ = numpy.array([numpy.log10(data_k), numpy.log10(data_Pk)])
			self._storeObject(table_name, table_)
			if self.print_info:
				print("Cosmology.matterPowerSpectrum: Lookup table completed.")	
			
			interpolator = self._getStoredObject(table_name, interpolator = True, inverse = inverse)

		return interpolator

	###############################################################################################

	def matterPowerSpectrum(self, k, Pk_source = 'eh98', derivative = False):
		"""
		The matter power spectrum at a scale k.
		
		By default, the power spectrum is computed using the transfer function approximation of 
		Eisenstein & Hu 1998 (``eh98``). Alternatively, the user can choose their version 
		without BAO (``eh98smooth``), or a user-submitted file with an arbitrary name. In that case,
		a file with two columns (k, P(k)) must be placed in the storage directory, and be named
		``matterpower_<cosmology_name>_<Pk_source>``, e.g. ``matterpower_planck13_camb``.
		
		The Eisenstein & Hu 1998 approximation is accurate to about 1%, and the interpolation 
		introduces errors significantly smaller than that.
		
		Parameters
		-------------------------------------------------------------------------------------------
		k: array_like
			The wavenumber k (in comoving h/Mpc), where :math:`10^{-20} < k < 10^{20}`; can be a 
			number or a numpy array.
		Pk_source: str
			Either ``eh98``, ``eh98smooth``, or the name of a user-defined table.
		derivative: bool
			If False, return P(k). If True, return :math:`d \log(P) / d \log(k)`.
			
		Returns
		-------------------------------------------------------------------------------------------
		Pk: array_like
			The matter power spectrum (or its logarithmic derivative if ``derivative == True``); has 
			the same dimensions as k.

		See also
		-------------------------------------------------------------------------------------------
		transferFunctionEH98: The transfer function according to Eisenstein & Hu 1998.
		transferFunctionEH98Smooth: The transfer function according to Eisenstein & Hu 1998, without the BAO features.
		"""
			
		if self.interpolation and (Pk_source == 'eh98' or Pk_source == 'eh98smooth'):
			
			# Load lookup-table
			interpolator = self._matterPowerSpectrumInterpolator(Pk_source)
			
			# If the requested radius is outside the range, give a detailed error message.
			k_req = numpy.min(k)
			if k_req < self.k_Pk[0]:
				msg = "k = %.2e is too small (min. k = %.2e)" % (k_req, self.k_min_Pk)
				raise Exception(msg)
		
			k_req = numpy.max(k)
			if k_req > self.k_Pk[-1]:
				msg = "k = %.2e is too large (max. k = %.2e)" % (k_req, self.k_max_Pk)
				raise Exception(msg)

			if derivative:
				Pk = interpolator(numpy.log10(k), nu = 1)
			else:
				Pk = interpolator(numpy.log10(k))
				Pk = 10**Pk
			
		else:
			
			if derivative > 0:
				raise Exception("Derivative can only be evaluated if interpolation == True.")

			if Utilities.isArray(k):
				Pk = k * 0.0
				for i in range(len(k)):
					Pk[i] = self._matterPowerSpectrumExact(k[i], Pk_source = Pk_source, ignore_norm = False)
			else:
				Pk = self._matterPowerSpectrumExact(k, Pk_source = Pk_source, ignore_norm = False)

		return Pk
	
	###############################################################################################

	def filterFunction(self, filt, k, R, no_oscillation = False):
		"""
		The filter function for the variance in Fourier space. 
		
		This function is dimensionless, the input units are k in comoving h/Mpc and R in comoving 
		Mpc/h. Please see the documentation of the :func:`sigma` function for details.

		Parameters
		-------------------------------------------------------------------------------------------
		filt: str
			Either ``tophat`` or ``gaussian``.
		k: float
			A wavenumber k (in comoving h/Mpc).
		R: float
			A radius R (in comoving Mpc/h).
		no_oscillation: bool
			If True, omit the sine / cosine part of the tophat filter; for internal use only.
			
		Returns
		-------------------------------------------------------------------------------------------
		filter: float
			The value of the filter function.
		"""
		
		x = k * R
		
		if filt == 'tophat':
			
			if no_oscillation:
				if x < 1.0:
					ret = 1.0
				else:
					ret = x**-2
			else:
				if x < 1E-3:
					ret = 1.0
				else:
					ret = 3.0 / x**3 * (numpy.sin(x) - x * numpy.cos(x))
				
		elif filt == 'gaussian':
			ret = numpy.exp(-x**2)
		
		else:
			msg = "Invalid filter, %s." % (filt)
			raise Exception(msg)
			
		return ret

	###############################################################################################

	def _sigmaExact(self, R, j = 0, filt = 'tophat', Pk_source = 'eh98', exact_Pk = False, ignore_norm = False):

		# -----------------------------------------------------------------------------------------
		def logIntegrand(lnk, Pk_interpolator, test = False):
			
			k = numpy.exp(lnk)
			W = self.filterFunction(filt, k, R, no_oscillation = test)
			
			if exact_Pk or (not self.interpolation):
				Pk = self._matterPowerSpectrumExact(k, Pk_source = Pk_source, ignore_norm = ignore_norm)
			else:
				Pk = 10**Pk_interpolator(numpy.log10(k))
			
			# One factor of k is due to the integration in log-k space
			ret = Pk * W**2 * k**3
			
			# Higher moment terms
			if j > 0:
				ret *= k**(2 * j)
			
			return ret

		# -----------------------------------------------------------------------------------------
		if filt == 'tophat' and j > 0:
			msg = "Higher-order moments of sigma are not well-defined for " + "tophat filter. Choose filt = 'gaussian' instead."
			raise Exception(msg)
	
		# For power-law cosmologies, we can evaluate sigma analytically. The exact expression 
		# has a dependence on n that in turn depends on the filter used, but the dependence 
		# on radius is simple and independent of the filter. Thus, we use sigma8 to normalize
		# sigma directly. 
		if self.power_law:
			
			n = self.power_law_n + 2 * j
			if n <= -3.0:
				msg = "n + 2j must be > -3 for the variance to converge in a power-law cosmology."
				raise Exception(msg)
			sigma2 = R**(-3 - n) / (8.0**(-3 - n) / self.sigma8**2)
			sigma = numpy.sqrt(sigma2)
			
		else:
			
			# If we are getting P(k) from a look-up table, it is a little more efficient to 
			# get the interpolator object and use it directly, rather than using the P(k) function.
			Pk_interpolator = None
			if (not exact_Pk) and self.interpolation:
				Pk_interpolator = self._matterPowerSpectrumInterpolator(Pk_source)
			
			# The infinite integral over k often causes trouble when the tophat filter is used. Thus,
			# we determine sensible limits and integrate over a finite volume. For tabled power 
			# spectra, we need to be careful not to exceed their limits.
			test_integrand_min = 1E-6
			test_k_min, test_k_max = self._matterPowerSpectrumLimits(Pk_source)
			test_k_min = max(test_k_min * 1.0001, 1E-7)
			test_k_max = min(test_k_max * 0.9999, 1E15)
			test_k = numpy.arange(numpy.log(test_k_min), numpy.log(test_k_max), 2.0)
			n_test = len(test_k)
			test_k_integrand = test_k * 0.0
			for i in range(n_test):
				test_k_integrand[i] = logIntegrand(test_k[i], Pk_interpolator)
			integrand_max = numpy.max(test_k_integrand)
			
			min_index = 0
			while test_k_integrand[min_index] < integrand_max * test_integrand_min:
				min_index += 1
				if min_index > n_test - 2:
					msg = "Could not find lower integration limit."
					raise Exception(msg)

			min_index -= 1
			max_index = min_index + 1
			while test_k_integrand[max_index] > integrand_max * test_integrand_min:
				max_index += 1	
				if max_index == n_test:
					msg = "Could not find upper integration limit."
					raise Exception(msg)
	
			args = Pk_interpolator
			sigma2, _ = scipy.integrate.quad(logIntegrand, test_k[min_index], test_k[max_index], \
						args = args, epsabs = 0.0, epsrel = self.accuracy_sigma, limit = 100)
			sigma = numpy.sqrt(sigma2 / 2.0 / math.pi**2)
		
		if numpy.isnan(sigma):
			msg = "Result is nan (cosmology %s, filter %s, R %.2e, j %d." % (self.name, filt, R, j)
			raise Exception(msg)
			
		return sigma
	
	###############################################################################################

	# Return a spline interpolator for sigma(R) or R(sigma) if inverse == True. Generally, sigma(R) 
	# should be evaluated using the sigma() function below, but for some performance-critical 
	# operations it is faster to obtain the interpolator directly from this function.If the lookup-
	# table does not exist yet, create it. For sigma, we use a very particular binning scheme. At 
	# low R, sigma is a very smooth function, and very wellapproximated by a spline interpolation 
	# between few points. Around the BAO scale, we need a higher resolution. Thus, the bins are 
	# assigned in reverse log(log) space.

	def _sigmaInterpolator(self, j, Pk_source, filt, inverse):
		
		table_name = 'sigma%d_%s_%s_%s' % (j, self.name, Pk_source, filt)
		interpolator = self._getStoredObject(table_name, interpolator = True, inverse = inverse)
		
		if interpolator is None:
			if self.print_info:
				print("Cosmology.sigma: Computing lookup table.")
			max_log = numpy.log10(self.R_max_sigma)
			log_range = max_log - numpy.log10(self.R_min_sigma)
			max_loglog = numpy.log10(log_range + 1.0)
			loglog_width = max_loglog / self.R_Nbins_sigma
			R_loglog = numpy.arange(0.0, max_loglog + loglog_width, loglog_width)
			log_R = max_log - 10**R_loglog[::-1] + 1.0
			data_R = 10**log_R
			data_sigma = data_R * 0.0
			for i in range(len(data_R)):
				data_sigma[i] = self._sigmaExact(data_R[i], j = j, filt = filt, Pk_source = Pk_source)
			table_ = numpy.array([numpy.log10(data_R), numpy.log10(data_sigma)])
			self._storeObject(table_name, table_)
			if self.print_info:
				print("Cosmology.sigma: Lookup table completed.")

			interpolator = self._getStoredObject(table_name, interpolator = True, inverse = inverse)
	
		return interpolator

	###############################################################################################
	
	def sigma(self, R, j = 0, z = 0.0, inverse = False, derivative = False, Pk_source = 'eh98', filt = 'tophat'):
		"""
		The rms variance of the linear density field on a scale R, :math:`\\sigma(R)`.
		
		The variance and its higher moments are defined as the integral
		
		.. math::
			\\sigma^2(R,z) = \\frac{1}{2 \\pi^2} \\int_0^{\\infty} k^2 k^{2j} P(k,z) |\\tilde{W}(kR)|^2 dk

		where :math:`\\tilde{W}(kR)$` is the Fourier transform of the :func:`filterFunction`, and 
		:math:`P(k,z) = D_+^2(z)P(k,0)` is the :func:`matterPowerSpectrum`. By default, the power 
		spectrum is computed using the transfer function approximation of Eisenstein & Hu 1998 
		(``eh98``) which is accurate to about 1%. The integration and interpolation introduce errors 
		smaller than that. Higher moments of the variance (such as :math:`\sigma_1`, 
		:math:`\sigma_2` etc) can be computed by setting j > 0 (see Bardeen et al. 1986). For the
		higher moments, the interpolation error increases to up to ~0.5%.
		Furthermore, the logarithmic derivative of :math:`\sigma(R)` can be evaluated by setting 
		``derivative == True``.
		
		Parameters
		-------------------------------------------------------------------------------------------
		R: array_like
			The radius of the filter in comoving Mpc/h, where :math:`10^{-12} < R < 10^3`; can be 
			a number or a numpy array.
		j: integer
			The order of the integral. j = 0 corresponds to the variance, j = 1 to the same integral 
			with an extra :math:`k^2` term etc; see Bardeen et al. 1986 for mathematical details.
		z: float
			Redshift; for z > 0, :math:`\sigma(R)` is multiplied by the linear growth factor.
		filt: str
			Either ``tophat`` or ``gaussian``. Higher moments (j > 0) can only be computed for the 
			gaussian filter.
		Pk_source: str
			Either ``eh98``, ``eh98smooth``, or the name of a user-supplied table.
		inverse: bool
			If True, compute :math:`R(\sigma)` rather than :math:`\sigma(R)`. For internal use.
		derivative: bool
			If True, return the logarithmic derivative, :math:`d \log(\sigma) / d \log(R)`, or its
			inverse, :math:`d \log(R) / d \log(\sigma)` if ``inverse == True``.
			
		Returns
		-------------------------------------------------------------------------------------------
		sigma: array_like
			The rms variance; has the same dimensions as R. If inverse and/or derivative are True, 
			the inverse, derivative, or derivative of the inverse are returned. If j > 0, those 
			refer to higher moments.

		See also
		-------------------------------------------------------------------------------------------
		matterPowerSpectrum: The matter power spectrum at a scale k.
		"""
		
		if self.interpolation:
			interpolator = self._sigmaInterpolator(j, Pk_source, filt, inverse)
			
			if not inverse:
	
				# If the requested radius is outside the range, give a detailed error message.
				R_req = numpy.min(R)
				if R_req < self.R_min_sigma:
					M_min = 4.0 / 3.0 * math.pi * self.R_min_sigma**3 * self.rho_m(0.0) * 1E9
					msg = "R = %.2e is too small (min. R = %.2e, min. M = %.2e)" \
						% (R_req, self.R_min_sigma, M_min)
					raise Exception(msg)
			
				R_req = numpy.max(R)
				if R_req > self.R_max_sigma:
					M_max = 4.0 / 3.0 * math.pi * self.R_max_sigma**3 * self.rho_m(0.0) * 1E9
					msg = "R = %.2e is too large (max. R = %.2e, max. M = %.2e)" \
						% (R_req, self.R_max_sigma, M_max)
					raise Exception(msg)
	
				if derivative:
					ret = interpolator(numpy.log10(R), nu = 1)
				else:
					ret = interpolator(numpy.log10(R))
					ret = 10**ret
					if z > 1E-5:
						ret *= self.growthFactor(z)

			else:
				
				sigma_ = R
				if z > 1E-5:
					sigma_ /= self.growthFactor(z)

				# Get the limits in sigma from storage, or compute and store them. Using the 
				# storage mechanism seems like overkill, but these numbers should be erased if 
				# the cosmology changes and sigma is re-computed.
				sigma_min = self._getStoredObject('sigma_min')
				sigma_max = self._getStoredObject('sigma_max')
				if sigma_min is None or sigma_min is None:
					knots = interpolator.get_knots()
					sigma_min = 10**numpy.min(knots)
					sigma_max = 10**numpy.max(knots)
					self._storeObject('sigma_min', sigma_min, persistent = False)
					self._storeObject('sigma_max', sigma_max, persistent = False)
				
				# If the requested sigma is outside the range, give a detailed error message.
				sigma_req = numpy.max(sigma_)
				if sigma_req > sigma_max:
					msg = "sigma = %.2e is too large (max. sigma = %.2e)" % (sigma_req, sigma_max)
					raise Exception(msg)
					
				sigma_req = numpy.min(sigma_)
				if sigma_req < sigma_min:
					msg = "sigma = %.2e is too small (min. sigma = %.2e)" % (sigma_req, sigma_min)
					raise Exception(msg)
				
				# Interpolate to get R(sigma)
				if derivative: 
					ret = interpolator(numpy.log10(sigma_), nu = 1)					
				else:
					ret = interpolator(numpy.log10(sigma_))
					ret = 10**ret

		else:
			
			if inverse:
				raise Exception('R(sigma), and thus massFromPeakHeight(), cannot be evaluated with interpolation == False.')
			if derivative:
				raise Exception('Derivative of sigma cannot be evaluated if interpolation == False.')

			if Utilities.isArray(R):
				ret = R * 0.0
				for i in range(len(R)):
					ret[i] = self._sigmaExact(R[i], j = j, filt = filt, Pk_source = Pk_source)
			else:
				ret = self._sigmaExact(R, j = j, filt = filt, Pk_source = Pk_source)
			if z > 1E-5:
				ret *= self.growthFactor(z)
		
		return ret

	###############################################################################################
	
	def peakHeight(self, M, z, filt = 'tophat', Pk_source = 'eh98', deltac_const = True):
		"""
		Peak height, :math:`\\nu`, given a halo mass.
		
		Peak height is defined as :math:`\\nu \equiv \delta_c / \sigma(M)`. See the documentation 
		of the :func:`sigma` function for details on filters, power spectrum source etc.
		
		Parameters
		-------------------------------------------------------------------------------------------
		M: array_like
			Halo mass in :math:`M_{\odot}/h`; can be a number or a numpy array.
		z: float
			Redshift.
		filt: str
			Either ``tophat`` or ``gaussian``.
		Pk_source: str
			Either ``eh98``, ``eh98smooth``, or the name of a user-supplied table.
		deltac_const: bool
			If True, the function returns the constant top-hat model collapse overdensity. If False,
			a correction due to the ellipticity of halos is applied.

		Returns
		-------------------------------------------------------------------------------------------
		nu: array_like
			Peak height; has the same dimensions as M.

		See also
		-------------------------------------------------------------------------------------------
		massFromPeakHeight: Halo mass from peak height, :math:`\\nu`.
		sigma: The rms variance of the linear density field on a scale R, :math:`\sigma(R)`.
		"""
				
		R = self.lagrangianR(M)
		sigma = self.sigma(R, z = z, filt = filt, Pk_source = Pk_source)
		nu = self.collapseOverdensity(deltac_const, sigma) / sigma

		return nu
	
	###############################################################################################

	def massFromPeakHeight(self, nu, z, filt = 'tophat', Pk_source = 'eh98', deltac_const = True):
		"""
		Halo mass from peak height, :math:`\\nu`.
		
		Peak height is defined as :math:`\\nu \equiv \delta_c / \sigma(M)`. See the documentation 
		of the :func:`sigma` function for details on filters, power spectrum source etc.
		
		Parameters
		-------------------------------------------------------------------------------------------
		nu: array_like
			Peak height; can be a number or a numpy array.
		z: float
			Redshift.
		filt: str
			Either ``tophat`` or ``gaussian``.
		Pk_source: str
			Either ``eh98``, ``eh98smooth``, or the name of a user-supplied table.
		deltac_const: bool
			If True, the function returns the constant top-hat model collapse overdensity. If False,
			a correction due to the ellipticity of halos is applied.

		Returns
		-------------------------------------------------------------------------------------------
		M: array_like
			Mass in :math:`M_{\odot}/h`; has the same dimensions as nu.

		See also
		-------------------------------------------------------------------------------------------
		peakHeight: Peak height, :math:`\\nu`, given a halo mass.
		sigma: The rms variance of the linear density field on a scale R, :math:`\sigma(R)`.
		"""
		
		sigma = self.collapseOverdensity(deltac_const = deltac_const) / nu
		R = self.sigma(sigma, z = z, filt = filt, Pk_source = Pk_source, inverse = True)
		M = self.lagrangianM(R)
		
		return M
	
	###############################################################################################
	
	def nonLinearMass(self, z, filt = 'tophat', Pk_source = 'eh98'):
		"""
		The non-linear mass, :math:`M^*`.
		
		:math:`M^*` is the mass for which the variance is equal to the collapse threshold, i.e.
		:math:`\sigma(M^*) = \delta_c` and thus :math:`\\nu(M^*) = 1`. See the documentation 
		of the :func:`sigma` function for details on filters, power spectrum source etc.
		
		Parameters
		-------------------------------------------------------------------------------------------
		z: float
			Redshift.
		filt: str
			Either ``tophat`` or ``gaussian``.
		Pk_source: str
			Either ``eh98``, ``eh98smooth``, or the name of a user-supplied table.

		Returns
		-------------------------------------------------------------------------------------------
		Mstar: float
			The non-linear mass in :math:`M_{\odot}/h`.

		See also
		-------------------------------------------------------------------------------------------
		massFromPeakHeight: Halo mass from peak height, :math:`\\nu`.
		sigma: The rms variance of the linear density field on a scale R, :math:`\sigma(R)`.
		"""
				
		return self.massFromPeakHeight(1.0, z = z, filt = filt, Pk_source = Pk_source, deltac_const = True)

	###############################################################################################
	# Peak curvature routines
	###############################################################################################
	
	# Get the mean peak curvature, <x>, at fixed nu from the integral of Bardeen et al. 1986 
	# (BBKS). Note that this function is approximated very well by the _peakCurvatureApprox() 
	# function below.
	
	def _peakCurvatureExact(self, nu, gamma):
	
		# Equation A15 in BBKS. 
		
		def curvature_fx(x):
	
			f1 = math.sqrt(5.0 / 2.0) * x
			t1 = scipy.special.erf(f1) + scipy.special.erf(f1 / 2.0)
	
			b0 = math.sqrt(2.0 / 5.0 / math.pi)
			b1 = 31.0 * x ** 2 / 4.0 + 8.0 / 5.0
			b2 = x ** 2 / 2.0 - 8.0 / 5.0
			t2 = b0 * (b1 * math.exp(-5.0 * x ** 2 / 8.0) + b2 * math.exp(-5.0 * x ** 2 / 2.0))
	
			res = (x ** 3 - 3.0 * x) * t1 / 2.0 + t2
	
			return res
	
		# Equation A14 in BBKS, minus the normalization which is irrelevant here. If we need the 
		# normalization, the Rstar parameter also needs to be passed.
		
		def curvature_Npk(x, nu, gamma):
	
			#norm = math.exp(-nu**2 / 2.0) / (2 * math.pi)**2 / Rstar**3
			norm = 1.0
			fx = curvature_fx(x)
			xstar = gamma * nu
			g2 = 1.0 - gamma ** 2
			exponent = -(x - xstar) ** 2 / (2.0 * g2)
			res = norm * fx * math.exp(exponent) / math.sqrt(2.0 * math.pi * g2)
	
			return res
	
		# Average over Npk
		
		def curvature_Npk_x(x, nu, gamma):
			return curvature_Npk(x, nu, gamma) * x
	
		args = nu, gamma
		norm, _ = scipy.integrate.quad(curvature_Npk, 0.0, numpy.infty, args, epsrel = 1E-10)
		integ, _ = scipy.integrate.quad(curvature_Npk_x, 0.0, numpy.infty, args, epsrel = 1E-10)
		xav = integ / norm
	
		return xav
	
	###############################################################################################
	
	# Wrapper for the function above which takes tables of sigmas. This form can be more convenient 
	# when computing many different nu's. 
	
	def _peakCurvatureExactFromSigma(self, sigma0, sigma1, sigma2, deltac_const = True):
	
		nu = self.collapseOverdensity(deltac_const, sigma0) / sigma0
		gamma = sigma1 ** 2 / sigma0 / sigma2
	
		x = nu * 0.0
		for i in range(len(nu)):
			x[i] = self._peakCurvatureExact(nu[i], gamma[i])
	
		return nu, gamma, x
	
	###############################################################################################
	
	# Get peak curvature from the approximate formula in BBKS. This approx. is excellent over the 
	# relevant range of nu.
	
	def _peakCurvatureApprox(self, nu, gamma):
	
		# Compute theta according to Equation 6.14 in BBKS
		g = gamma
		gn = g * nu
		theta1 = 3.0 * (1.0 - g ** 2) + (1.216 - 0.9 * g ** 4) * numpy.exp(-g * gn * gn / 8.0)
		theta2 = numpy.sqrt(3.0 * (1.0 - g ** 2) + 0.45 + (gn / 2.0) ** 2) + gn / 2.0
		theta = theta1 / theta2
	
		# Equation 6.13 in BBKS
		x = gn + theta
		
		# Equation 6.15 in BBKS
		nu_tilde = nu - theta * g / (1.0 - g ** 2)
	
		return theta, x, nu_tilde
	
	###############################################################################################
	
	# Wrapper for the function above which takes tables of sigmas. This form can be more convenient 
	# when computing many different nu's. For convenience, various intermediate numbers are 
	# returned as well.
	
	def _peakCurvatureApproxFromSigma(self, sigma0, sigma1, sigma2, deltac_const = True):
	
		nu = self.collapseOverdensity(deltac_const, sigma0) / sigma0
		gamma = sigma1**2 / sigma0 / sigma2
		
		theta, x, nu_tilde = self._peakCurvatureApprox(nu, gamma)
		
		return nu, gamma, x, theta, nu_tilde
	
	###############################################################################################
	
	def peakCurvature(self, M, z, filt = 'gaussian', Pk_source = 'eh98', \
					deltac_const = True, exact = False):
		"""
		The average curvature of peaks for a halo mass M.
		
		In a Gaussian random field, :math:`\delta`, the peak height is defined as 
		:math:`\delta / \\sigma` where :math:`\\sigma = \\sigma_0` is the rms variance. The 
		curvature of the field is defined as :math:`x = -\\nabla^2 \delta / \\sigma_2` where 
		:math:`\\sigma_2` is the second moment of the variance.
		
		This function computes the average curvature of peaks in a Gaussian random field, <x>,
		according to Bardeen et al. 1986 (BBKS), for halos of a certain mass M. This mass is 
		converted to a Lagrangian scale R, and thus the variance and its moments. The evaluation
		can be performed by integration of Equation A14 in BBKS (if ``exact == True``), or using their
		fitting function in Equation 6.13 (if ``exact == False``). The fitting function is excellent over
		the relevant range of peak heights. 
		
		Parameters
		-------------------------------------------------------------------------------------------
		M: array_like
			Mass in in :math:`M_{\odot}/h`; can be a number or a numpy array.
		z: float
			Redshift.
		filt: str
			Either ``tophat`` or ``gaussian``.
		Pk_source: str
			Either ``eh98``, ``eh98smooth``, or the name of a user-supplied table.
		deltac_const: bool
			If ``True``, the function returns the constant top-hat model collapse overdensity. If 
			``False``, a correction due to the ellipticity of halos is applied.
		exact: bool
			If ``True``, evaluate the integral exactly; if ``False``, use the BBKS approximation.	

		Returns
		-------------------------------------------------------------------------------------------
		nu: array_like
			Peak height; has the same dimensions as M.
		gamma: array_like
			An intermediate parameter, :math:`\\gamma = \\sigma_1^2 / (\\sigma_0 \\sigma_2)` (see
			Equation 4.6a in BBKS); has the same dimensions as M.
		x: array_like
			The mean peak curvature for halos of mass M (note the caveat discussed above); has the 
			same dimensions as M.
		theta: array_like
			An intermediate parameter (see Equation 6.14 in BBKS; only returned if ``exact == False``); 
			has the same dimensions as M.
		nu_tilde: array_like
			The modified peak height (see Equation 6.15 in BBKS; only returned if ``exact == False``); 
			has the same dimensions as M.
		
		Warnings
		-------------------------------------------------------------------------------------------		
		While peak height quantifies how high a fluctuation over the background a halo is, peak
		curvature tells us something about the shape of the initial peak. However, note the 
		cloud-in-cloud problem (BBKS): not all peaks end up forming halos, particularly small
		peaks will often get swallowed by other peaks. Thus, the average peak curvature is not
		necessarily equal to the average curvature of peaks that form halos.		
		"""

		R = self.lagrangianR(M)
		sigma0 = self.sigma(R, j = 0, z = z, filt = filt, Pk_source = Pk_source)
		sigma1 = self.sigma(R, j = 1, z = z, filt = filt, Pk_source = Pk_source)
		sigma2 = self.sigma(R, j = 2, z = z, filt = filt, Pk_source = Pk_source)
	
		if exact:
			return self._peakCurvatureExactFromSigma(sigma0, sigma1, sigma2, deltac_const = deltac_const)
		else:
			return self._peakCurvatureApproxFromSigma(sigma0, sigma1, sigma2, deltac_const = deltac_const)

	###############################################################################################

	def _correlationFunctionExact(self, R, Pk_source = 'eh98'):
		
		f_cut = 0.001

		# -----------------------------------------------------------------------------------------
		# The integrand is exponentially cut off at a scale 1000 * R.
		def integrand(k, R, Pk_source, Pk_interpolator):
			
			if self.interpolation:
				Pk = 10**Pk_interpolator(numpy.log10(k))
			else:
				Pk = self.matterPowerSpectrumExact(k, Pk_source)

			ret = Pk * k / R * numpy.exp(-(k * R * f_cut)**2)
			
			return ret

		# -----------------------------------------------------------------------------------------
		# If we are getting P(k) from a look-up table, it is a little more efficient to 
		# get the interpolator object and use it directly, rather than using the P(k) function.
		Pk_interpolator = None
		if self.interpolation:
			Pk_interpolator = self._matterPowerSpectrumInterpolator(Pk_source)

		# Use a Clenshaw-Curtis integration, i.e. an integral weighted by sin(kR). 
		k_min = 1E-6 / R
		k_max = 10.0 / f_cut / R
		args = R, Pk_source, Pk_interpolator
		xi, _ = scipy.integrate.quad(integrand, k_min, k_max, args = args, epsabs = 0.0, \
					epsrel = self.accuracy_xi, limit = 100, weight = 'sin', wvar = R)
		xi /= 2.0 * math.pi**2

		if numpy.isnan(xi):
			msg = 'Result is nan (cosmology %s, R %.2e).' % (self.name, R)
			raise Exception(msg)

		return xi
	
	###############################################################################################

	# Return a spline interpolator for the correlation function, xi(R). Generally, xi(R) should be 
	# evaluated using the correlationFunction() function below, but for some performance-critical 
	# operations it is faster to obtain the interpolator directly from this function.

	def _correlationFunctionInterpolator(self, Pk_source):

		table_name = 'correlation_%s_%s' % (self.name, Pk_source)
		interpolator = self._getStoredObject(table_name, interpolator = True)
		
		if interpolator is None:
			if self.print_info:
				print("correlationFunction: Computing lookup table. This may take a few minutes, please do not interrupt.")
			
			data_R = numpy.zeros((numpy.sum(self.R_xi_Nbins) + 1), numpy.float)
			n_regions = len(self.R_xi_Nbins)
			k_computed = 0
			for i in range(n_regions):
				log_min = numpy.log10(self.R_xi[i])
				log_max = numpy.log10(self.R_xi[i + 1])
				log_range = log_max - log_min
				bin_width = log_range / self.R_xi_Nbins[i]
				if i == n_regions - 1:
					data_R[k_computed:k_computed + self.R_xi_Nbins[i] + 1] = \
						10**numpy.arange(log_min, log_max + bin_width, bin_width)
				else:
					data_R[k_computed:k_computed + self.R_xi_Nbins[i]] = \
						10**numpy.arange(log_min, log_max, bin_width)
				k_computed += self.R_xi_Nbins[i]
			
			data_xi = data_R * 0.0
			for i in range(len(data_R)):
				data_xi[i] = self._correlationFunctionExact(data_R[i], Pk_source = Pk_source)
			table_ = numpy.array([data_R, data_xi])
			self._storeObject(table_name, table_)
			if self.print_info:
				print("correlationFunction: Lookup table completed.")
			interpolator = self._getStoredObject(table_name, interpolator = True)
		
		return interpolator

	###############################################################################################

	def correlationFunction(self, R, z = 0.0, derivative = False, Pk_source = 'eh98'):
		"""
		The linear matter-matter correlation function at radius R.
		
		The linear correlation function is defined as 
		
		.. math::
			\\xi(R) = \\frac{1}{2 \\pi^2} \\int_0^\\infty k^2 P(k) \\frac{\\sin(kR)}{kR} dk
		
		where P(k) is the :func:`matterPowerSpectrum`. The integration, as well as the 
		interpolation routine, are accurate to ~1-2% over the range :math:`10^{-3} < R < 500`. 
		
		Parameters
		-------------------------------------------------------------------------------------------
		R: array_like
			The radius of the filter in comoving Mpc/h; can be a number or a numpy array.
		z: float
			Redshift
		derivative: bool
			If ``derivative == True``, the linear derivative :math:`d \\xi / d R` is returned.
		Pk_source: str
			Either ``eh98``, ``eh98smooth``, or the name of a user-supplied table.

		Returns
		-------------------------------------------------------------------------------------------
		xi: array_like
			The correlation function, or its derivative; has the same dimensions as R.

		See also
		-------------------------------------------------------------------------------------------
		matterPowerSpectrum: The matter power spectrum at a scale k.
		"""
		
		if self.interpolation:
			
			# Load lookup-table
			interpolator = self._correlationFunctionInterpolator(Pk_source)
				
			# If the requested radius is outside the range, give a detailed error message.
			R_req = numpy.min(R)
			if R_req < self.R_xi[0]:
				msg = 'R = %.2e is too small (min. R = %.2e)' % (R_req, self.R_xi[0])
				raise Exception(msg)
		
			R_req = numpy.max(R)
			if R_req > self.R_xi[-1]:
				msg = 'R = %.2e is too large (max. R = %.2e)' % (R_req, self.R_xi[-1])
				raise Exception(msg)
	
			# Interpolate to get xi(R). Note that the interpolation is performed in linear 
			# space, since xi can be negative.
			if derivative:
				ret = interpolator(R, nu = 1)
			else:
				ret = interpolator(R)
			
		else:

			if derivative:
				raise Exception('Derivative of xi cannot be evaluated if interpolation == False.')

			if Utilities.isArray(R):
				ret = R * 0.0
				for i in range(len(R)):
					ret[i] = self.correlationFunctionExact(R[i], Pk_source = Pk_source)
			else:
				ret = self.correlationFunctionExact(R, Pk_source = Pk_source)

		if not derivative and z > 1E-5:
			ret *= self.growthFactor(z)**2
		
		return	ret

###################################################################################################
# Setter / getter functions for cosmologies
###################################################################################################

def setCosmology(cosmo_name, params = None):
	"""
	Set a cosmology.
	
	This function provides a convenient way to create a cosmology object without setting the 
	parameters of the Cosmology class manually. See the Basic Usage section for examples.
	Whichever way the cosmology is set, the global variable is updated so that the :func:`getCurrent` 
	function returns the set cosmology.

	Parameters
	-----------------------------------------------------------------------------------------------
	cosmo_name: str
		The name of the cosmology.
	params: dictionary
		The parameters of the constructor of the Cosmology class.

	Returns
	-----------------------------------------------------------------------------------------------
	cosmo: Cosmology
		The created cosmology object.
	"""
	
	if 'powerlaw_' in cosmo_name:
		n = float(cosmo_name.split('_')[1])
		param_dict = cosmologies['powerlaw']
		param_dict['power_law'] = True
		param_dict['power_law_n'] = n
		if params is not None:
			param_dict = dict(param_dict.items() + params.items())
			
	elif cosmo_name in cosmologies:		
		param_dict = cosmologies[cosmo_name]
		if params is not None:
			param_dict = dict(param_dict.items() + params.items())
			
	else:
		if params is not None:
			param_dict = params.copy()
		else:
			msg = "Invalid cosmology (%s)." % (cosmo_name)
			raise Exception(msg)
		
	param_dict['name'] = cosmo_name
	cosmo = Cosmology(**(param_dict))
	setCurrent(cosmo)
	
	return cosmo

###################################################################################################

def addCosmology(cosmo_name, params):
	"""
	Add a set of cosmological parameters to the global list.
	
	After this function is executed, the new cosmology can be set using :func:`setCosmology` from 
	anywhere in the code.
	
	Parameters
	-----------------------------------------------------------------------------------------------
	cosmo_name: str
		The name of the cosmology.
	params: dictionary
		A set of parameters for the constructor of the Cosmology class.
	"""
	
	cosmologies[cosmo_name] = params
	
	return 

###################################################################################################

def setCurrent(cosmo):
	"""
	Set the current global cosmology to a cosmology object.
	
	Unlike :func:`setCosmology`, this function does not create a new cosmology object, but allows 
	the user to set a cosmology object to be the current cosmology. This can be useful when switching
	between cosmologies, since many routines use the :func:`getCurrent` routine to obtain the current
	cosmology.
	
	Parameters
	-----------------------------------------------------------------------------------------------
	cosmo: Cosmology
		The cosmology object to be set as the global current cosmology.
	"""

	global current_cosmo
	current_cosmo = cosmo
	
	return

###################################################################################################

def getCurrent():
	"""
	Get the current global cosmology.
	
	This function should be used whenever access to the cosmology is needed. By using the globally
	set cosmology, there is no need to pass cosmology objects around the code. If no cosmology is
	set, this function raises an Exception that reminds the user to set a cosmology.

	Returns
	-----------------------------------------------------------------------------------------------
	cosmo: Cosmology
		The current globally set cosmology. 
	"""
	
	if current_cosmo is None:
		raise Exception('Cosmology is not set.')

	return current_cosmo

###################################################################################################
