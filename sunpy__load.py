#!/usr/bin/env python
""" Useful loading routines for opening and doing basic processing of the SUNRISE fits files

Most/all loading tools here open the SUNRISE fits files to grab one or two basic fields, and 
return it to the use for further processing.  These routines are made use of in the sunpy__plot
routines.  The routines provided here are all very simple.  Many additions/extensions can
be included for further specific functionality.

For usage examples, see the plotting routines in sunpy__plot.
"""
import numpy as np
import os
import sys
import astropy.io.fits as fits
import cosmocalc
import pyfits 
import scipy as sp
import scipy.ndimage
import congrid
import matplotlib.pyplot as plt
import sunpy.sunpy__synthetic_image


__author__ = "Paul Torrey and Greg Snyder"
__copyright__ = "Copyright 2014, The Authors"
__credits__ = ["Paul Torrey", "Greg Snyder"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Paul Torrey"
__email__ = "ptorrey@cfa.harvard.edu"
__status__ = "Production"
if __name__ == '__main__':    #code to execute if called from command-line
    pass    #do nothing 



def my_fits_open(filename):
    if (not os.path.exists(filename)):
        print "file not found:", filename
        sys.exit()
    return fits.open(filename)


def load_broadband_image(filename,band=0,camera=0):
  """ Loads an idealized sunrise broadband image for a specified fits file, band, and camera.
      The band can be specified as a number or a string (must match the "band_names")		"""

  band_images = load_all_broadband_images(filename,camera=camera)
  band_names  = load_broadband_names(filename)

  if type(band) is int:
    return_image = band_images[band,:,:]
  else:
    band = (((band_names == band).nonzero())[0])[0]
    return_image = band_images[band,:,:]

  return return_image


def load_fov(filename):
    hdulist = my_fits_open(filename)
    data = hdulist['CAMERA0-PARAMETERS'].header['linear_fov']
    hdulist.close()
    return data


def load_broadband_names(filename):
    hdulist = my_fits_open(filename)
    name_array = hdulist['FILTERS'].data.field(0)
    hdulist.close()
    return name_array


def load_broadband_fast_names(filename):
    print " "
    print "WARNING: fast NAMES HAVE BEEN HARD-CODED; CHECK OUTPUT BELOW FOR CONSISTENCY!!!"

    sunrise_names = load_broadband_names(filename)
    fast_names    = sunrise_names 			## initial guess

    fast_names[2] = "SDSS/u.dat"
    fast_names[3] = "SDSS/g.dat"
    fast_names[4] = "SDSS/r.dat"
    fast_names[5] = "SDSS/i.dat"
    fast_names[6] = "SDSS/z.dat"

    for index in range(len(fast_names)):
        print sunrise_names[index]
        print fast_names[index]
        print " "

    return fast_names


def load_broadband_effective_wavelengths(filename,band=None):
  if (not os.path.exists(filename)):
    print "file not found:", filename
    sys.exit()

  hdulist = fits.open(filename)
  name_array = hdulist['FILTERS'].data['lambda_eff']
  hdulist.close()
  if band != None:
    if type(band) is int:
      name_array = name_array[band]
    else:
      band_names  = load_broadband_names(filename)
      band_index = (((band_names == band).nonzero())[0])[0]
      name_array = name_array[band_index]

  band=None
  return name_array


def load_all_broadband_images(filename,camera=0):
  if (not os.path.exists(filename)):
    print "file not found:", filename
    sys.exit()

  camera_string = 'CAMERA'+str(camera)+'-BROADBAND-NONSCATTER'
  hdulist = fits.open(filename)
  data = hdulist[camera_string].data

  data[ data < 1e-20 ] = 1e-20

  hdulist.close()
  return data


def load_broadband_image(filename,band=0,camera=0):
  band_images = load_all_broadband_images(filename,camera=camera)
  band_names  = load_broadband_names(filename)
 
  if type(band) is int:
    return_image = band_images[band,:,:]
  else:
    band = (((band_names == band).nonzero())[0])[0]
    return_image = band_images[band,:,:]

  return return_image


def load_all_broadband_photometry(filename,camera=0):
  if (not os.path.exists(filename)):
    print "file not found:", filename
    return 0

  hdulist = fits.open(filename)
  data = hdulist['FILTERS'].data['AB_mag_nonscatter0']  
  return data
 
 
def load_all_broadband_apparent_magnitudes(filename,camera=0,dist=4e8):
    dist_modulus = 5.0 * ( np.log10(dist) - 1.0 )
    apparent_magnitudes = dist_modulus + load_all_broadband_photometry(filename,camera=camera)
    return apparent_magnitudes



def load_redshift(filename):
  if (not os.path.exists(filename)):
    print "file not found:", filename
    sys.exit()

  hdulist = fits.open(filename)
  redshift = hdulist[1].header['REDSHIFT']
  hdulist.close()
  return redshift


def load_sed_lambda(filename):
  hdulist = my_fits_open(filename)
  lambda_array = hdulist['INTEGRATED_QUANTITIES'].data['lambda  ']
  hdulist.close()
  return lambda_array


def load_sed_l_lambda(filename):
  hdulist = my_fits_open(filename)
  l_lambda_array = hdulist['INTEGRATED_QUANTITIES'].data['L_lambda']
  hdulist.close()
  return l_lambda_array



