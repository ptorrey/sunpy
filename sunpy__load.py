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
try:
    import astropy.io.fits as fits
    print "loaded astropy.io.fits"
except:
    try:
        import pyfits as fits
        print "loaded pyfits"
    except:
        print "Error: Unable to access PyFITS or AstroPy modules.\n\n"+"With root access, add PyFITS to your site-packages with:\n\n"+"% pip install pyfits\n"+"or\n"+"% easy_install pyfits\n\n"+"or download at: www.stsci.edu/institute/software_hardware/pyfits/Download\n"+"where additional installation options and instructions can be found."

#import astropy.io.fits as fits
import cosmocalc			# http://cxc.harvard.edu/contrib/cosmocalc/
import scipy as sp
import scipy.ndimage
import sunpy.sunpy__synthetic_image


__author__ = "Paul Torrey and Greg Snyder"
__copyright__ = "Copyright 2014, The Authors"
__credits__ = ["Paul Torrey", "Greg Snyder"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Paul Torrey"
__email__ = "ptorrey@cfa.harvard.edu"
__status__ = "Production"
if __name__ == '__main__':
    pass

speedoflight_m  = 2.99e8

def my_fits_open(filename):
    if (not os.path.exists(filename)):
        print "file not found:", filename
        sys.exit()
    return fits.open(filename)


def load_broadband_image(filename,band=0, **kwargs):
  """ Loads an idealized sunrise broadband image for a specified fits file, band, and camera.
      The band can be specified as a number or a string (must match the "band_names")		"""

  band_images = load_all_broadband_images(filename, **kwargs)
  band_names  = load_broadband_names(filename)

  if type(band) is int:
    return_image = band_images[band,:,:]
  else:
    band = (((band_names == band).nonzero())[0])[0]
    return_image = band_images[band,:,:]

  return return_image

def load_broadband_magnitude(filename, band=0, **kwargs):
  band_mags = load_all_broadband_photometry(filename, **kwargs) 
  band_names  = load_broadband_names(filename)
  if type(band) is not int:
      band = int( np.where([this_band == band for this_band in band_names])[0][0]  )
  return_val = band_mags[band]

  return return_val

def load_resolved_broadband_magnitudes(filename, band=0, camera=0, **kwargs):
    """ this is a little trickier b/c in W/m/m^2/str.  First convert to abs mag, then dist correction """
    band_names  = load_broadband_names(filename)
    if type(band) is not int:
        band = int( np.where([this_band == band for this_band in band_names])[0][0]  )

#    if type(band) is not int:
#      band = int( (((band_names == band).nonzero())[0])[0] )

    image = load_broadband_image(filename,band=band,camera=camera)	# in W/m/m^2/str  shape = [n_band, n_pix, n_pix]
    mag   = load_broadband_magnitude(filename, band=band, camera=camera)

    n_pixels = image.shape[1]

    hdulist = fits.open(filename)
    lambda_eff = hdulist['FILTERS'].data['lambda_eff']
    this_lambda = lambda_eff[band]
    
    to_nu                     = ((this_lambda**2 ) / (speedoflight_m)) #* pixel_area_in_str
    to_microjanskies          = (1.0e6) * to_nu * (1.0e26)                 # 1 muJy/str (1Jy = 1e-26 W/m^2/Hz)
    image *=  to_microjanskies              # to microjanskies / str

    pixel_in_kpc           = load_fov(filename)  / n_pixels
    pixel_in_sr = (1e3 * pixel_in_kpc / 10.0)**2
    image *=  pixel_in_sr                                      # in muJy
    image /= 1e6                                               # in Jy                         
    tot_img_in_Jy = np.sum(image)                               # total image flux in Jy
    abmag = -2.5 * np.log10(tot_img_in_Jy / 3631 )
    print abmag
    return abmag


def load_fov(filename):
    hdulist = my_fits_open(filename)
    data = hdulist['CAMERA0-PARAMETERS'].header['linear_fov']
    hdulist.close()
    return data

def load_camera_angles(filename,camera=0):
    hdulist = my_fits_open(filename)
    theta = hdulist['CAMERA'+str(camera)+'-PARAMETERS'].header['theta']
    phi   = hdulist['CAMERA'+str(camera)+'-PARAMETERS'].header['phi']
    hdulist.close()
    return theta,phi


def load_broadband_names(filename):
    hdulist = my_fits_open(filename)
    name_array = hdulist['FILTERS'].data.field(0)
    hdulist.close()
    name_array = [ s.replace(" ", "") for s in name_array ]
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


def load_all_broadband_images(filename,camera=0,openlist=None,use_nonscatter=True):
  if (not os.path.exists(filename)):
    print "file not found:", filename
    sys.exit()

  if use_nonscatter is True:
      camera_string = 'CAMERA'+str(camera)+'-BROADBAND-NONSCATTER'
  else:
      camera_string = 'CAMERA'+str(camera)+'-BROADBAND'
      
  if openlist is None:
      openlist = fits.open(filename,memmap=False)
      data = openlist[camera_string].data
      #openlist.close()
      print "### Sunpy: opening broadband list: ", openlist.filename()
  else:
      openfn = openlist.filename()
      assert openfn==filename
      data = openlist[camera_string].data

  data[ data < 1e-20 ] = 1e-20

  return data,openlist


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

  cst = str(camera)
  hdulist = fits.open(filename)
  data = hdulist['FILTERS'].data['AB_mag_nonscatter'+cst]  
  return data
 
 
def load_integrated_broadband_apparent_magnitudes(filename,camera=0,dist=4e8):
    """ this is fairly easy b/c already in abs mag.  Only need to do distance correction """
    dist_modulus = 5.0 * ( np.log10(dist) - 1.0 )
    apparent_magnitudes = dist_modulus + load_all_broadband_photometry(filename,camera=camera)
    return apparent_magnitudes


def load_resolved_broadband_apparent_magnitudes(filename, redshift, camera=0, **kwargs):
    """ this is a little trickier b/c in W/m/m^2/str.  First convert to abs mag, then dist correction """
    images = load_all_broadband_images(filename,camera=0)        # in W/m/m^2/str  shape = [n_band, n_pix, n_pix]
    mags   = load_all_broadband_photometry(filename, camera=0)

    n_pixels = images.shape[1]

    hdulist = fits.open(filename)
    lambda_eff = hdulist['FILTERS'].data['lambda_eff']
    for index,this_lambda in enumerate(lambda_eff):
        to_nu                     = ((this_lambda**2 ) / (speedoflight_m)) #* pixel_area_in_str
        to_microjanskies          = (1.0e6) * to_nu * (1.0e26)                 # 1 muJy/str (1Jy = 1e-26 W/m^2/Hz)
        images[index,:,:] = images[index,:,:] * to_microjanskies              # to microjanskies / str
    
    pixel_in_kpc           = load_fov(filename)  / n_pixels
    pixel_in_sr = (1e3 * pixel_in_kpc / 10.0)**2
    images *=  pixel_in_sr                 			# in muJy
    images /= 1e6						# in Jy				
    for index,this_lambda in enumerate(lambda_eff):
        tot_img_in_Jy = np.sum(images[index,:,:])           		  	# total image flux in Jy
        abmag = -2.5 * np.log10(tot_img_in_Jy / 3631 )
	if True:
            print "the ab magnitude of this image is :"+str(abmag)+"  "+str(mags[index])
	    print abmag/mags[index], abmag - mags[index]

        print index, np.sum(images[index,:,:])

    images = -2.5 * np.log10( images / 3631 )			# abmag in each pixel

    dist = (cosmocalc.cosmocalc(redshift, H0=70.4, WM=0.2726, WV=0.7274))['DL_Mpc'] * 1e6
    dist_modulus = 5.0 * ( np.log10(dist) - 1.0 )
    apparent_magnitudes = dist_modulus + images
    
    return apparent_magnitudes

#    apparent_magnitudes = dist_modulus + 



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

# these options only with for sunrise with rad. transfer.  Not for current Illustris images #
def load_sed_l_lambda_with_rt(filename, camera=0):
  hdulist = my_fits_open(filename)
  l_lambda_array = hdulist['INTEGRATED_QUANTITIES'].data['L_lambda_out'+str(camera)]
  hdulist.close()
  return l_lambda_array

def load_sed_l_lambda_scatter(filename, camera=0):
  hdulist = my_fits_open(filename)
  l_lambda_array = hdulist['INTEGRATED_QUANTITIES'].data['L_lambda_scatter'+str(camera)]
  hdulist.close()
  return l_lambda_array

def load_sed_l_lambda_nonscatter(filename, camera=0):
  hdulist = my_fits_open(filename)
  l_lambda_array = hdulist['INTEGRATED_QUANTITIES'].data['L_lambda_nonscatter'+str(camera)]
  hdulist.close()
  return l_lambda_array

def load_sed_l_lambda_ir(filename, camera=0):
  hdulist = my_fits_open(filename)
  l_lambda_array = hdulist['INTEGRATED_QUANTITIES'].data['L_lambda_ir'+str(camera)]
  hdulist.close()
  return l_lambda_array
#===============================================================================#


#===============================================================================#
def load_stellar_mass_map(filename,camera=0):
  hdulist = fits.open(filename)
  camera_string = 'CAMERA'+str(camera)+'-AUX'
  aux_image = hdulist[camera_string].data
  map = aux_image[4,:,:]
  return map

def load_mass_weighted_stellar_age_map(filename,camera=0):
  hdulist = fits.open(filename)
  camera_string = 'CAMERA'+str(camera)+'-AUX'
  aux_image = hdulist[camera_string].data
  map = aux_image[7,:,:]
  return map

def load_stellar_metal_map(filename,camera=0):
  hdulist = fits.open(filename)
  camera_string = 'CAMERA'+str(camera)+'-AUX'
  aux_image = hdulist[camera_string].data
  map = aux_image[5,:,:]
  return map
