"""
routines for plotting images from the SUNRISE data output


Example usage:

Dependencies:

"""


import numpy as np
import os
import sys
import sunpy__load			# used for noiseless images, for which we can return the image input directly
import sunpy__synthetic_image		# used for images with noise, pixel scaling, etc.
import matplotlib.pyplot as plt
import img_scale
import gc




__author__ = "Paul Torrey and Greg Snyder"
__copyright__ = "Copyright 2014, The Authors"
__credits__ = ["Paul Torrey", "Greg Snyder"]
__license__ = "GPL"
__version__ = "1.0"
__maintainer__ = "Paul Torrey"
__email__ = "ptorrey@mit.harvard.edu"
__status__ = "Production"
if __name__ == '__main__':    #code to execute if called from command-line
    pass    #do nothing 




def plot_synthetic_sdss_gri(filename,camera=0,savefile=None,scale_min=1e-3, lupton_alpha=1, lupton_Q=1e-10, psf_fwhm_arcsec=None,
                                        include_background=True ):
    """ routine for plotting synthetic sdss gri images from Illustris idealized images including appropriate pixel scaling, noise, etc.  """

    if (not os.path.exists(filename)):
        print "file not found:", filename
        sys.exit()

    if savefile==None:
        savefile = 'synthetic_image_'+filename[filename.index('broadband_')+10:filename.index('.fits')]+'.png'

    rp, img = return_synthetic_sdss_gri_img(filename, camera=camera, scale_min=scale_min, lupton_alpha=lupton_alpha, lupton_Q=lupton_Q, psf_fwhm_arcsec=psf_fwhm_arcsec,
                                        include_background=include_background)

    my_save_image(img, savefile, opt_text='rp='+str(rp))
    del img
    gc.collect()

def plot_sdss_gri(filename,camera=0,savefile='./sdss_gri.png',scale_min=0.1,scale_max=50, non_linear=0.5 ):
    """ routine for plotting synthetic sdss gri images from Illustris idealized images *without* additional image effects """

    if (not os.path.exists(filename)):
        print "file not found:", filename
        sys.exit()

    img = return_sdss_gri_img(filename,camera=camera,scale_min=scale_min,scale_max=scale_max, non_linear=non_linear)
    my_save_image(img, savefile)
    del img
    gc.collect()


def return_synthetic_sdss_gri_img(filename, camera=0, lupton_alpha=0.5, lupton_Q=0.5, scale_min=1e-4, psf_fwhm_arcsec=None, 
					include_background=True):
    seed=int(filename[filename.index('broadband_')+10:filename.index('.fits')])

    b_image, rp    = sunpy__synthetic_image.build_synthetic_image(filename, 'g_SDSS.res', r_petro_kpc=None, seed=(seed*(camera+1)), camera=camera, psf_fwhm_arcsec=psf_fwhm_arcsec, include_background=include_background)
    g_image, dummy = sunpy__synthetic_image.build_synthetic_image(filename, 'r_SDSS.res', r_petro_kpc=rp,   seed=(seed*(camera+1)), camera=camera, psf_fwhm_arcsec=psf_fwhm_arcsec, include_background=include_background)
    r_image, dummy = sunpy__synthetic_image.build_synthetic_image(filename, 'i_SDSS.res', r_petro_kpc=rp,   seed=(seed*(camera+1)), camera=camera, psf_fwhm_arcsec=psf_fwhm_arcsec, include_background=include_background)
    n_pixels = r_image.shape[0]
    img = np.zeros((n_pixels, n_pixels, 3), dtype=float)

    b_image *= 1.1
    g_image *= 1.0
    r_image *= 1.0
 
    I = (r_image + g_image + b_image)/3
    val = np.arcsinh( lupton_alpha * lupton_Q * (I - scale_min))/lupton_Q

    img[:,:,0] = r_image * val / I
    img[:,:,1] = g_image * val / I
    img[:,:,2] = b_image * val / I

    print "Lupton Q     "+str(lupton_Q)
    print "scale_min    "+str(scale_min)
    print "img min/max/mean "+str(img.min())+"  "+str(img.max())+"  "+str(img.mean())
    print " "

    maxrgbval = np.amax(img, axis=2)
 
    changeind = maxrgbval > 1.0
    img[changeind,0] = img[changeind,0]/maxrgbval[changeind]
    img[changeind,1] = img[changeind,1]/maxrgbval[changeind]
    img[changeind,2] = img[changeind,2]/maxrgbval[changeind]

    minrgbval = np.amin(img, axis=2)
    changeind = minrgbval < 0.0
    img[changeind,0] = 0
    img[changeind,1] = 0
    img[changeind,2] = 0

    changind = I < 0
    img[changind,0] = 0
    img[changind,1] = 0
    img[changind,2] = 0
    img[img<0] = 0

    print "img min/max/mean "+str(img.min())+"  "+str(img.max())+"  "+str(img.mean())

    del b_image, g_image, r_image, I, val
    gc.collect()

    img[img<0] = 0
    return rp, img

def return_sdss_gri_img(filename,camera=0,scale_min=0.1,scale_max=50,size_scale=1.0, non_linear=0.5):
    if (not os.path.exists(filename)):
        print "file not found:", filename
        sys.exit()

    b_image = sunpy__load.load_broadband_image(filename,band='g_SDSS.res',camera=camera) * 0.7
    g_image = sunpy__load.load_broadband_image(filename,band='r_SDSS.res',camera=camera) * 1.0
    r_image = sunpy__load.load_broadband_image(filename,band='i_SDSS.res',camera=camera) * 1.4
    n_pixels = r_image.shape[0]
    img = np.zeros((n_pixels, n_pixels, 3), dtype=float)

    img[:,:,0] = img_scale.asinh(r_image, scale_min=scale_min, scale_max=scale_max,non_linear=non_linear)
    img[:,:,1] = img_scale.asinh(g_image, scale_min=scale_min, scale_max=scale_max,non_linear=non_linear)
    img[:,:,2] = img_scale.asinh(b_image, scale_min=scale_min, scale_max=scale_max,non_linear=non_linear)
    img[img<0] = 0

    del b_image, g_image, r_image
    gc.collect()

    return img


def return_h_band_img(filename,camera=0,scale_min=0.1,scale_max=50,size_scale=1.0):
    if (not os.path.exists(filename)):
        print "file not found:", filename
        sys.exit()

    image = sunpy__load.load_broadband_image(filename,band='H_Johnson.res', camera=camera) 
    n_pixels = image.shape[0]
    img = np.zeros((n_pixels, n_pixels), dtype=float)
    img[:,:] = img_scale.asinh(image, scale_min=scale_min, scale_max=scale_max,non_linear=0.5)
    img[img<0] = 0
    return img


def return_johnson_uvk_img(filename,camera=0,scale_min=0.1,scale_max=50,size_scale=1.0):
    if (not os.path.exists(filename)):
        print "file not found:", filename
        sys.exit()

    b_effective_wavelength = sunpy__load.load_broadband_effective_wavelengths(filename,band="U_Johnson.res")
    g_effective_wavelength = sunpy__load.load_broadband_effective_wavelengths(filename,band="V_Johnson.res")
    r_effective_wavelength = sunpy__load.load_broadband_effective_wavelengths(filename,band="K_Johnson.res")

    b_image = sunpy__load.load_broadband_image(filename,band='U_Johnson.res',camera=camera) * b_effective_wavelength / g_effective_wavelength * 2.5
    g_image = sunpy__load.load_broadband_image(filename,band='V_Johnson.res',camera=camera) * g_effective_wavelength / g_effective_wavelength 
    r_image = sunpy__load.load_broadband_image(filename,band='K_Johnson.res',camera=camera) * r_effective_wavelength / g_effective_wavelength * 1.5 

    n_pixels
    img = np.zeros((n_pixels, n_pixels, 3), dtype=float)
    img[:,:,0] = img_scale.asinh(r_image, scale_min=scale_min, scale_max=scale_max,non_linear=0.5)
    img[:,:,1] = img_scale.asinh(g_image, scale_min=scale_min, scale_max=scale_max,non_linear=0.5)
    img[:,:,2] = img_scale.asinh(b_image, scale_min=scale_min, scale_max=scale_max,non_linear=0.5)
    img[img<0] = 0
    return img



def my_save_image(img, savefile, opt_text=None):
    if img.shape[0] >1:
        fig = plt.figure(figsize=(1,1))
        ax = fig.add_subplot(111)
        imgplot = ax.imshow(img,origin='lower')
        plt.axis('off')

        if not opt_text==None:
            ax.text(img.shape[0]/2.0, img.shape[0]/2.0, opt_text, ha='center',va='center', color='white', fontsize=4)

        fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
        fig.savefig(savefile, dpi=img.shape[0]-1)
        fig.clf()
        plt.close()
        del img
        gc.collect()
