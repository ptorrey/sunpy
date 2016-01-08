"""
routines for plotting images from the SUNRISE data output


Example usage:

import sunpy.sunpy__plot as sunplot

sunplot.plot_synthetic_sdss_gri('broadband_1234.fits')
sunplot.plot_sdss_gri('broadband_1234.fits')

Dependencies:
  numpy
  matplotlib
"""

import numpy as np
import os
import sys
import sunpy.sunpy__load as sunpy__load			# used for noiseless images, for which we can return the image input directly
import sunpy__synthetic_image		# used for images with noise, pixel scaling, etc.
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
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



def plot_spectrum(filelist, savefile='spectrum.pdf', fontsize=14, ymin=None, ymax=None, **kwargs):
    """ routine for plotting the spectra for a list of files """
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111)
    print filelist
    for file in filelist:
	print file
	wavelength = sunpy__load.load_sed_lambda(file)
        sed        = sunpy__load.load_sed_l_lambda(file)
        ax.plot(wavelength,wavelength * sed, label=file)

    ax.set_xlabel(r'$\lambda$ (m)',fontsize=fontsize)
    ax.set_ylabel(r'$\lambda$L${}_{\lambda}$ (W)', fontsize=fontsize)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(fontsize)

    ax.legend
    ax.set_xscale('log')
    ax.set_yscale('log')
    if (ymin != None) and (ymax != None):
	ax.set_ylim([ymin, ymax])

    fig.subplots_adjust(left=0.18, right=0.95, top=0.95, bottom=0.12, wspace=0.0, hspace=0.0)
    fig.savefig(savefile)





def plot_synthetic_sdss_gri(filename, savefile='syn_sdss_gri.png', **kwargs):
    """ routine for plotting synthetic sdss gri images from Illustris idealized images including appropriate pixel scaling, noise, etc.  """

    rp, img = return_synthetic_sdss_gri_img(filename, **kwargs)
    my_save_image(img, savefile)
    del img
    gc.collect()

def plot_synthetic_hst(filename, savefile='syn_hst.png', **kwargs):	#mass_string=None, full_fov=None, cut_bad_pixels=False, **kwargs):
    """ routine for plotting synthetic hst gri images from Illustris idealized images including appropriate pixel scaling, noise, etc.  """

    rp, img = return_synthetic_hst_img(filename, **kwargs)	#return_synthetic_sdss_gri_img(filename, **kwargs)
    my_save_image(img, savefile, **kwargs)	#top_opt_text=mass_string, full_fov=full_fov, cut_bad_pixels=cut_bad_pixels, **kwargs)
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
				seed_boost=1.0,
			 	r_petro_kpc=None,
				**kwargs):

    fail_flag=True		# looks for "bad" backgrounds, and tells us to try again
    n_iter = 1
    while(fail_flag and (n_iter < 2)):
        fail_flag=False
	try:
            seed=int(filename[filename.index('broadband_')+10:filename.index('.fits')])*(n_iter)*seed_boost
        except:
	    try:
                seed=int(filename[filename.index('.fits')-3:filename.index('.fits')])*(n_iter)*seed_boost
	    except:
		seed=1234

        n_iter+=1


        b_image, rp, the_used_seed,this_fail_flag    = sunpy__synthetic_image.build_synthetic_image(filename, 'g_SDSS.res', 
				seed=seed,
				r_petro_kpc=r_petro_kpc, 
				fix_seed=False,
				**kwargs)
        if(this_fail_flag):
	  fail_flag=True

        g_image, dummy, the_used_seed,this_fail_flag = sunpy__synthetic_image.build_synthetic_image(filename, 'r_SDSS.res', 
				seed=the_used_seed,
				r_petro_kpc=rp,
				fix_seed=True, 
				**kwargs)
        if(this_fail_flag):
          fail_flag=True

        r_image, dummy, the_used_seed, this_fail_flag = sunpy__synthetic_image.build_synthetic_image(filename, 'i_SDSS.res', 
                                seed=the_used_seed,
				r_petro_kpc=rp,
				fix_seed=True, 
				**kwargs)
        if(this_fail_flag):
            fail_flag=True

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


    del b_image, g_image, r_image, I, val
    gc.collect()

    img[img<0] = 0
    return rp, img


def return_synthetic_hst_img(filename,
                                lupton_alpha=0.5, lupton_Q=0.5, scale_min=1e-4,
                                b_fac=1.0, g_fac=1.0, r_fac=1.0, max=1.0, dynrng=1e3,
                                **kwargs):

    try:
       seed=int(filename[filename.index('broadband_')+10:filename.index('.fits')])
    except:
       try:
          seed=int(filename[filename.index('broadband_red_')+14:filename.index('.fits')])
       except:
          seed=int(filename[filename.index('broadband_rest_')+15:filename.index('.fits')])
    

    b_image, rp, dummy, dummy    = sunpy__synthetic_image.build_synthetic_image(filename, 22,		#25,
                                seed=seed, fix_seed=True,
                                #r_petro_kpc=None,
                                **kwargs)
    dummy, dummy, dummy, dummy  = sunpy__synthetic_image.build_synthetic_image(filename, 25,
                                seed=seed, fix_seed=True,
                                #r_petro_kpc=rp,
                                **kwargs)

    g_image, dummy, dummy, dummy = sunpy__synthetic_image.build_synthetic_image(filename, 26,
                                seed=seed, fix_seed=True,
                                #r_petro_kpc=rp,
                                **kwargs)
    r_image, dummy, dummy, dummy = sunpy__synthetic_image.build_synthetic_image(filename, 27,
                                seed=seed, fix_seed=True,
                                #r_petro_kpc=rp,
                                **kwargs)

    n_pixels = r_image.shape[0]
    img = np.zeros((n_pixels, n_pixels, 3), dtype=float)

    b_image *= b_fac
    g_image *= g_fac
    r_image *= r_fac

    if 0:
       I = (r_image + g_image + b_image)/3
       val = np.arcsinh( lupton_alpha * lupton_Q * (I - scale_min))/lupton_Q
       I[ I < 1e-8 ] = 1e20               # from below, this effectively sets the pixel to 0

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

    else:
        img[:,:,0] = np.log10(r_image )  # / max=1.0, dynrng=1e3,
	img[:,:,1] = np.log10(g_image )
	img[:,:,2] = np.log10(b_image )

        img -= np.log10( max/dynrng ) 
        img /= np.log10( dynrng ) 

        img[img>1] = 1
	I = img
	val = img

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

    img[:,:,0] = asinh(r_image, scale_min=scale_min, scale_max=scale_max,non_linear=non_linear)
    img[:,:,1] = asinh(g_image, scale_min=scale_min, scale_max=scale_max,non_linear=non_linear)
    img[:,:,2] = asinh(b_image, scale_min=scale_min, scale_max=scale_max,non_linear=non_linear)
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
    img[:,:] = asinh(image, scale_min=scale_min, scale_max=scale_max,non_linear=0.5)
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
    img[:,:,0] = asinh(r_image, scale_min=scale_min, scale_max=scale_max,non_linear=0.5)
    img[:,:,1] = asinh(g_image, scale_min=scale_min, scale_max=scale_max,non_linear=0.5)
    img[:,:,2] = asinh(b_image, scale_min=scale_min, scale_max=scale_max,non_linear=0.5)
    img[img<0] = 0
    return img


def return_stellar_mass_img(filename, camera=0, scale_min=1e8, scale_max=1e10, size_scale=1.0, non_linear=1e8):
    image = sunpy__load.load_stellar_mass_map(filename, camera=camera)
    n_pixels = image.shape[0]
    img = np.zeros((n_pixels, n_pixels), dtype=float)
    img[:,:] = asinh(image, scale_min=scale_min, scale_max=scale_max, non_linear=non_linear)
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

def my_save_image(img, savefile, opt_text=None, top_opt_text=None, full_fov=None, cut_bad_pixels=False, zoom=None, save_indiv_bands=True, **kwargs):
    if img.shape[0] >1:
	n_pixels_save = img.shape[0]


	if zoom!=None:
	    n_pixels_current = img.shape[0]
	    center = np.floor(img.shape[0]/(2.0))
	    n_pixels_from_center  = np.floor(img.shape[0]/(2.0*zoom))
	    img = img[center - n_pixels_from_center:center+n_pixels_from_center, center - n_pixels_from_center:center+n_pixels_from_center]
	    print "resized image from "+str(n_pixels_current)+" to "+str(img.shape[0])+" pixels"
	    if full_fov!=None:
		full_fov /= (1.0*zoom)

	if cut_bad_pixels:
    	    cut_index=-1
	    print np.mean(img[cut_index:, cut_index:])
	    while( np.mean(img[cut_index:, cut_index:]) == 0 ):
	        print cut_index
	        cut_index -= 1
	    img = img[:cut_index,:cut_index]


        fig = plt.figure(figsize=(1,1))
        ax = fig.add_subplot(111)
        imgplot = ax.imshow(img,origin='lower', interpolation='nearest')
        plt.axis('off')

        if not opt_text==None:
            ax.text(img.shape[0]/2.0, img.shape[0]/2.0, opt_text, ha='center',va='center', color='white', fontsize=4)

        if not top_opt_text==None:
            ax.text(img.shape[0]/2.0, 0.9*img.shape[0], top_opt_text, ha='center',va='center', color='white', fontsize=4)

	if not full_fov==None:
	    bar_size_in_kpc    = np.round(full_fov/5.0)
	    pixel_size_in_kpc  = full_fov / (1.0 * img.shape[0])
	    bar_size_in_pixels = bar_size_in_kpc / pixel_size_in_kpc
	    center = img.shape[0]/2.0
	    
	    ax.text( center, 0.15 * img.shape[0], str(bar_size_in_kpc)+" kpc", color='w', fontsize=4, ha='center', va='center')
	    ax.plot( [center-bar_size_in_pixels/2.0, center+bar_size_in_pixels/2.0] , [0.1*img.shape[0], 0.1*img.shape[0]], lw=2,color='w')
	    ax.set_xlim([0,img.shape[0]-1])
	    ax.set_ylim([0,img.shape[0]-1])


        fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
        fig.savefig(savefile, dpi=n_pixels_save)
        fig.clf()
        plt.close()


	if save_indiv_bands:
	    print img.shape
	    for iii in np.arange(3):
		fig = plt.figure(figsize=(1,1))
                ax = fig.add_subplot(111)
		print img[:,:,iii].min(), img[:,:,iii].max()

                imgplot = ax.imshow(img[:,:,iii],origin='lower', interpolation='nearest', cmap = cm.Greys, vmin=0, vmax=1)
                plt.axis('off')
                fig.subplots_adjust(left=0.0, right=1.0, top=1.0, bottom=0.0)
                fig.savefig(savefile+"_band_"+str(iii)+'.png', dpi=n_pixels_save)
                fig.clf()
                plt.close()
        del img
        gc.collect()





def asinh(inputArray, scale_min=None, scale_max=None, non_linear=2.0):
        imageData=np.array(inputArray, copy=True)

        if scale_min == None:
                scale_min = imageData.min()
        if scale_max == None:
                scale_max = imageData.max()
        factor = np.arcsinh((scale_max - scale_min)/non_linear)
        indices0 = np.where(imageData < scale_min)
        indices1 = np.where((imageData >= scale_min) & (imageData <= scale_max))
        indices2 = np.where(imageData > scale_max)
        imageData[indices0] = 0.0
        imageData[indices2] = 1.0
        imageData[indices1] = np.arcsinh((imageData[indices1] - scale_min)/non_linear)/factor

        return imageData
