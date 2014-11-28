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




def plot_synthetic_sdss_gri(filename, savefile='syn_sdss_gri.png', **kwargs):
    """ routine for plotting synthetic sdss gri images from Illustris idealized images including appropriate pixel scaling, noise, etc.  """

    rp, img = return_synthetic_sdss_gri_img(filename, **kwargs)
    my_save_image(img, savefile)
    del img
    gc.collect()


def plot_sdss_gri(filename, savefile='./sdss_gri.png', **kwargs):
    """ routine for plotting synthetic sdss gri images from Illustris idealized images *without* additional image effects """

    img = return_sdss_gri_img(filename, **kwargs)
    my_save_image(img, savefile)
    del img
    gc.collect()


def return_synthetic_sdss_gri_img(filename, 
				lupton_alpha=0.5, lupton_Q=0.5, scale_min=1e-4, 
                                b_fac=0.7, g_fac=1.0, r_fac=1.3,
				**kwargs):

    seed=int(filename[filename.index('broadband_')+10:filename.index('.fits')])
    b_image, rp    = sunpy__synthetic_image.build_synthetic_image(filename, 'g_SDSS.res', 
				seed=seed,
				r_petro_kpc=None, 
				**kwargs)
    g_image, dummy = sunpy__synthetic_image.build_synthetic_image(filename, 'r_SDSS.res', 
				seed=seed,
				r_petro_kpc=rp, 
				**kwargs)
    r_image, dummy = sunpy__synthetic_image.build_synthetic_image(filename, 'i_SDSS.res', 
                                seed=seed,
				r_petro_kpc=rp, 
				**kwargs)

    n_pixels = r_image.shape[0]
    img = np.zeros((n_pixels, n_pixels, 3), dtype=float)
 
    b_image *= b_fac
    g_image *= g_fac
    r_image *= r_fac
 
    I = (r_image + g_image + b_image)/3
    val = np.arcsinh( lupton_alpha * lupton_Q * (I - scale_min))/lupton_Q
    I[ I < 1e-6 ] = 1e100		# from below, this effectively sets the pixel to 0

    img[:,:,0] = r_image * val / I
    img[:,:,1] = g_image * val / I
    img[:,:,2] = b_image * val / I

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
    print " "

    del b_image, g_image, r_image, I, val
    gc.collect()

    img[img<0] = 0
    return rp, img

def return_synthetic_hst_img(filename,
                                lupton_alpha=0.5, lupton_Q=0.5, scale_min=1e-4,
                                b_fac=1.0, g_fac=1.0, r_fac=1.0,
                                **kwargs):

    seed=int(filename[filename.index('broadband_')+10:filename.index('.fits')])
    b_image, rp    = sunpy__synthetic_image.build_synthetic_image(filename, 22,		#25,
                                seed=seed,
                                r_petro_kpc=None,
                                **kwargs)
    g_image, dummy = sunpy__synthetic_image.build_synthetic_image(filename, 26,
                                seed=seed,
                                r_petro_kpc=rp,
                                **kwargs)
    r_image, dummy = sunpy__synthetic_image.build_synthetic_image(filename, 27,
                                seed=seed,
                                r_petro_kpc=rp,
                                **kwargs)

    print b_image.min(), b_image.max()
    print g_image.min(), g_image.max()
    print r_image.min(), r_image.max()


    n_pixels = r_image.shape[0]
    img = np.zeros((n_pixels, n_pixels, 3), dtype=float)

    b_image *= b_fac
    g_image *= g_fac
    r_image *= r_fac

    I = (r_image + g_image + b_image)/3
    val = np.arcsinh( lupton_alpha * lupton_Q * (I - scale_min))/lupton_Q
    I[ I < 1e-8 ] = 1e20               # from below, this effectively sets the pixel to 0

    print val.min(), val.max()
    print I.min(), I.max()

    img[:,:,0] = r_image * val / I
    img[:,:,1] = g_image * val / I
    img[:,:,2] = b_image * val / I

    print " "
    print " "
    print img.min(), img.max()

    print " "
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
    print " "

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


def return_stellar_mass_img(filename, camera=0, scale_min=1e8, scale_max=1e10, size_scale=1.0, non_linear=1e8):
    image = sunpy__load.load_stellar_mass_map(filename, camera=camera)
    n_pixels = image.shape[0]
    img = np.zeros((n_pixels, n_pixels), dtype=float)
    img[:,:] = img_scale.asinh(image, scale_min=scale_min, scale_max=scale_max, non_linear=non_linear)
    img[img<0] = 0
    return img

def return_mass_weighted_age_img(filename, camera=0, scale_min=None, scale_max=None, size_scale=1.0):
    image = sunpy__load.load_mass_weighted_stellar_age_map(filename, camera=camera)
    image += 1e5
    return image

def return_stellar_metal_img(filename, camera=0, scale_min=None, scale_max=None, size_scale=1.0, non_linear=None):
    image1 = sunpy__load.load_stellar_mass_map(filename, camera=camera)
    image2 = sunpy__load.load_stellar_metal_map(filename, camera=camera)
    image = image2/image1
    image[image<0]      = 0     #image.min()
    image[image*0 != 0] = 0     #image.min()
    return image

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
