###################################################################################################
#
# Utilities.py 		        (c) Benedikt Diemer
#							University of Chicago
#     				    	bdiemer@oddjob.uchicago.edu
#
###################################################################################################

"""
Common routines for Colossus modules.
"""

import os
import numpy
import sys

###################################################################################################

def printLine():
	"""
	Print a line to the console.
	"""

	print('-------------------------------------------------------------------------------------')

	return

###################################################################################################

def getCacheDir(module = None):
	"""
	Get a directory for the persistent caching of data. The function attempts to locate the home 
	directory and (if necessary) create a .colossus sub-directory. In the rare case where that 
	fails, the location of this code file is used as a base directory.

	Parameters
	---------------------------
	module: string
		The name of the module that is requesting this cache directory. Each module has its own
		directory in order to avoid name conflicts.
	
	Returns
	-------
	path : string
		The cache directory.
	"""
	
	base_dir = getHomeDir()
	if base_dir is None:
		base_dir = getCodeDir()
	
	cache_dir = base_dir + '/.colossus/cache/'
	
	if module is not None:
		cache_dir += module + '/'

	if not os.path.exists(cache_dir):
		os.makedirs(cache_dir)
	
	return cache_dir

###################################################################################################

def getHomeDir():
	""" 
	Finds the home directory on this system.

	Returns
	-------
	path : string
		The home directory, or None if home cannot be found.
	"""
	
	def decodePath(path):
		return path.decode(sys.getfilesystemencoding())
	
	# There are basically two options for the operating system, either it's POSIX compatible of 
	# windows. POSIX includes UNIX, LINUX, Mac OS etc. The following choices were inspired by the 
	# astropy routine for finding a home directory.
	
	if os.name == 'posix':
		
		if 'HOME' in os.environ:
			home_dir = decodePath(os.environ['HOME'])
		else:
			raise Warning('Could not find HOME variable on POSIX-compatible operating system.')
			home_dir = None
	
	elif os.name == 'nt':
	
		if 'MSYSTEM' in os.environ and os.environ.get('HOME'):
			home_dir = decodePath(os.environ['HOME'])
		elif 'HOMESHARE' in os.environ:
			home_dir = decodePath(os.environ['HOMESHARE'])
		elif 'HOMEDRIVE' in os.environ and 'HOMEPATH' in os.environ:
			home_dir = os.path.join(os.environ['HOMEDRIVE'], os.environ['HOMEPATH'])
			home_dir = decodePath(home_dir)
		elif 'USERPROFILE' in os.environ:
			home_dir = decodePath(os.path.join(os.environ['USERPROFILE']))
		elif 'HOME' in os.environ:
			home_dir = decodePath(os.environ['HOME'])
		else:
			raise Warning('Could not find HOME directory on Windows system.')
			home_dir = None

	else:
	
		msg = 'Unknown operating system type, %s. Cannot find home directory.' % os.name
		raise Warning(msg)
		home_dir = None
	
	return home_dir

###################################################################################################

def getCodeDir():
	"""
	Returns the path to this code file.
	
	Returns
	-------
	path : string
		The code directory.
	"""
	
	path = os.path.dirname(os.path.realpath(__file__))

	return path

###################################################################################################

def isArray(var):
	"""
	Tests whether a variable var is iterable or not.

	Parameters
	---------------------------
	var: array_like
		Variable to be tested.
	
	Returns
	-------
	is_array: boolean
		Whether var is a numpy array or not.
	"""
	
	try:
		dummy = iter(var)
	except TypeError:
		is_array = False
	else:
		is_array = True
		
	return is_array

###################################################################################################

def getArray(var):
	"""
	Convert a variable to a numpy array, whether it already is one or not.

	Parameters
	-----------------------------------------------------------------------------------------------
	var: array_like
		Variable to be converted to an array.
	
	Returns
	-----------------------------------------------------------------------------------------------
	var_ret: numpy array
		A numpy array with one or more entries.
	is_array: boolean
		Whether var is a numpy array or not.
	"""
		
	is_array = isArray(var)
	if is_array:
		var_ret = var
	else:
		var_ret = numpy.array([var])
		
	return var_ret, is_array 

###################################################################################################
