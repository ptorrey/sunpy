#!/usr/bin/env python
""" Defines the synthetic_image class to add image realism to the idealized sunrise images.

The synthetic_image class defines a number of routines to take the original image and 
convolve it with some psf function, add sky noise, rebin to an appropate pixel scale 
(based on telescope), scale to an approparte image size (based on a petrosian radius 
calculation, and add background image (SDSS only supported at the moment).


The majority of the code in this file was developed by Greg Snyder and can be found in
	 Snyder et al., (2015), http://arxiv.org/abs/1502.07747 


"""
import numpy as np
import os
import sys
import math
import gc

try:
    import astropy.io.fits as fits
    print "loaded astropy.io.fits"
except:
    try:
        import pyfits as fits
        print "loaded pyfits"
    except:
        print "Error: Unable to access PyFITS or AstroPy modules."
        print "Add PyFITS to your site-packages with:"
        print "% pip install pyfits\n"
        print "  or  "
        print "% easy_install pyfits\n"
        print "  or  "
        print "download at: www.stsci.edu/institute/software_hardware/pyfits/Download\n"



import cosmocalc
import scipy as sp
import scipy.ndimage
import scipy.signal
import scipy.interpolate



try:
    import astropy.convolution.convolve as convolve ; CONVOLVE_TYPE='astropy'
    print "loaded astropy.convolution.convolve"
    from astropy.convolution import *
except:
    try:
        from scipy.signal import convolve2d as convolve ; CONVOLVE_TYPE='scipy'
        print "loaded scipy.signal.convolve2d; note that the astropy.convolution.convolve() function is preferred. There may be unexpected sub-pixel or off-by-one behavior with this scipy function."
    except:
        print "Error: Unable to access SciPy or AstroPy convolution modules."


import sunpy.sunpy__load
import time

import wget
import warnings

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

verbose=False


abs_dist        = 0.01         
erg_per_joule   = 1e7
speedoflight_m  = 2.99e8
m2_to_cm2       = 1.0e-4
n_arcsec_per_str = 4.255e10             # (radian per arc second)^2
n_pixels_galaxy_zoo = 424 


###########################################################
# SDSS background images created by Greg Snyder on 6/18/14#
# SDSS background obtained from:  data.sdss3.org/mosaics  #
# Ra = 175.0
# Dec = 30.0
# Size (deg) = 0.5
# Pixel Scale = 0.24 "/pixel
#
# HST backgrounds provided by Erica Nelson and Pascal Oesch 
# and integrated here by P. Torrey
# #########################################################

dl_base="http://www.illustris-project.org/files/backgrounds"
bg_base='./data/'
backgrounds = [	[], [], 		# GALEX 0 1
		[bg_base+'/SDSS_backgrounds/J113959.99+300000.0-u.fits'], 	# 2 SDSS-u 
		[bg_base+'/SDSS_backgrounds/J113959.99+300000.0-g.fits'], 	# 3 SDSS-g
		[bg_base+'/SDSS_backgrounds/J113959.99+300000.0-r.fits'], 	# 4 SDSS-r
		[bg_base+'/SDSS_backgrounds/J113959.99+300000.0-i.fits'], 	# 5 SDSS-i
		[bg_base+'/SDSS_backgrounds/J113959.99+300000.0-z.fits'], 	# 6 SDSS-z
		[], [], [], [],				# 7-8-9-10 IRAC
		[], [], [], [], [], [], [], [], [], [], 	# 11-12-13-14-15-16-17-18 JOHNSON/COUSINS + 2 mass
		[bg_base+'/HST_backgrounds/xdf_noise_F775W_30mas.fits'], #21	ACS-435
		[bg_base+'/HST_backgrounds/GOODSN_F606W.fits'], #22	ACS-606
		[bg_base+'/HST_backgrounds/xdf_noise_F775W_30mas.fits'], #23 	ACS-775
		[bg_base+'/HST_backgrounds/xdf_noise_F775W_30mas.fits'], #24	ACS-850
		[bg_base+'/HST_backgrounds/GOODSN_F125W.fits'],		 #25	f105w
	        [bg_base+'/HST_backgrounds/GOODSN_F125W.fits'],		 #26    f125w
		[bg_base+'/HST_backgrounds/GOODSN_F160W.fits'],		 #27	f160w
		[], [], [], [], [], [], [], []		# NIRCAM
		]

bg_zpt = {  "u_SDSS.res":[22.5],
            "g_SDSS.res":[22.5],
            "r_SDSS.res":[22.5],
            "i_SDSS.res":[22.5],
            "z_SDSS.res":[22.5]}



  #bg_zpt = [ [], [],                 # GALEX
#                        [22.5],
#                        [22.5],
#                        [22.5],
#                        [22.5],
#                        [22.5],
#                [], [], [], [],                         # 7-8-9-10 IRAC
#                [], [], [], [], [], [], [], [], [], [],         # 11-12-13-14-15-16-17-18 JOHNSON/COUSINS + 2 mass
#                [25.69],
#                [25.69],
#                [25.69],
#                [25.69],
#                [25.69],
#                [25.69],
#                [25.69],
#                [], [], [], [], [], [], [], []          # NIRCAM
#                ]


def build_synthetic_image(filename, band, r_petro_kpc=None, openlist=None, **kwargs):
    """ build a synthetic image from a SUNRISE fits file and return the image to the user """
    obj     	 = synthetic_image(filename, band=band, r_petro_kpc=r_petro_kpc, openlist=openlist, **kwargs)

    image = obj.bg_image.return_image()
    rp    = obj.r_petro_kpc
    seed  = obj.seed
    failed= obj.bg_failed
    fitsfn = obj.fitsfn
    openlist = obj.openlist

    del obj
    gc.collect() 
    
    return image, rp, seed, failed, fitsfn, openlist 

def load_resolved_broadband_apparent_magnitudes(filename, redshift, camera=0, seed=12345, n_bands=36, **kwargs):
    """ loads n_band x n_pix x n_pix image array with apparent mags for synthetic images """
    mags   = sunpy.sunpy__load.load_all_broadband_photometry(filename, camera=0)

    for band in np.arange(n_bands):
        obj          = synthetic_image(filename, band=int(band), seed=seed, redshift=redshift, **kwargs)    
	img = obj.bg_image.return_image()		#  muJy / str
	if band==0:
	    n_pixels = img.shape[0]
	    all_images = np.zeros( (n_bands, n_pixels, n_pixels ) )
	all_images[band, :, :] = img			# muJy / str

        pixel_in_sr = (1e3*obj.bg_image.pixel_in_kpc /10.0)**2
    
    all_images *=  pixel_in_sr / 1e6    	# in Jy

    for band in np.arange(n_bands):
        tot_img_in_Jy = np.sum(all_images[band,:,:])    # total image flux in Jy
        abmag = -2.5 * np.log10(tot_img_in_Jy / 3631 )
        if verbose:
            print "the ab magnitude of band "+str(band)+" is :"+str(abmag)+"  "+str(mags[band])
            print abmag/mags[band], abmag - mags[band]
            print " "

    all_images = -2.5 * np.log10( all_images / 3631 )                   # abmag in each pixel
    dist = (cosmocalc.cosmocalc(redshift, H0=70.4, WM=0.2726, WV=0.7274))['DL_Mpc'] * 1e6
    dist_modulus = 5.0 * ( np.log10(dist) - 1.0 )
    apparent_magnitudes = dist_modulus + all_images

    del mags, obj, img, n_pixels, all_images, pixel_in_sr, tot_img_in_Jy, abmag, dist, dist_modulus
    gc.collect()

    return apparent_magnitudes


class synthetic_image:
    """ main class for loading and manipulating SUNRISE data into real data format  """
    def __init__(self, 
            filename, band=0, camera=0,
            redshift=0.05,
            psf_fwhm_arcsec=1.0, pixelsize_arcsec=0.24,
            r_petro_kpc=None, save_fits=False,
            seed=None,
            add_background=True,
            add_psf=True,
            add_noise=True,
            rebin_phys=True,
            rebin_gz=False,
            n_target_pixels=n_pixels_galaxy_zoo,
            resize_rp=True,
            sn_limit=25.0,
            sky_sig=None,
            verbose=False,
            fix_seed=True,
	    bg_tag=None,
            bb_label='broadband_',
            output_label='',
            psf_fits=None,
            psf_pixsize_arcsec=None,
            psf_truncate_pixels=None,
            psf_hdu_num = 0,
            custom_fitsfile=None,
            bb_header=None,
            openlist=None,
            use_nonscatter=True,
            **kwargs):

        if (not os.path.exists(filename)):
            print "file not found:", filename
            sys.exit()

        start_time = time.time()
        self.filename  = filename
        self.cosmology = cosmology(redshift)
        self.telescope = telescope(psf_fwhm_arcsec, pixelsize_arcsec, psf_fits, psf_pixsize_arcsec, rebin_phys, add_psf, psf_truncate_pixels,psf_hdu_num)

        band_names  = sunpy.sunpy__load.load_broadband_names(filename)
        hdulist = fits.open(filename)
	
        if type(band) is not int:
	    band = int( np.where([this_band == band for this_band in band_names])[0][0]  )

        self.camera           = camera
        self.band             = band
        self.band_name        = band_names[band]
        if use_nonscatter is True:
            self.image_header     = hdulist['CAMERA'+str(camera)+'-BROADBAND-NONSCATTER'].header
        else:
            self.image_header     = hdulist['CAMERA'+str(camera)+'-BROADBAND'].header

        bb_header = self.image_header
        self.broadband_header = hdulist['BROADBAND'].header
        self.param_header     = hdulist['CAMERA'+str(camera)+'-PARAMETERS'].header
        self.int_quant_data   = hdulist['INTEGRATED_QUANTITIES'].data
        self.filter_data      = hdulist['FILTERS'].data
        self.lambda_eff       = (self.filter_data['lambda_eff'])[band]
        self.ewidth_lambda    = (self.filter_data['ewidth_lambda'])[band]
        self.ewidth_nu        = (self.filter_data['ewidth_nu'])[band]
        if use_nonscatter is True:
            self.sunrise_absolute_mag = (self.filter_data['AB_mag_nonscatter'+str(self.camera)])[band]
        else:
            self.sunrise_absolute_mag = (self.filter_data['AB_mag'+str(self.camera)])[band]
            
        hdulist.close()

#============= DECLARE ALL IMAGES HERE =================#
        self.sunrise_image  = single_image()		# orig sunrise image
        self.psf_image      = single_image()		# supersampled image + psf convolution
        self.rebinned_image = single_image()		# rebinned by appropriate pixel scale
        self.noisy_image    = single_image()		# noise added via gaussian draw
        self.nmag_image     = single_image()		# converted to nanomaggies units
        self.rp_image       = single_image()		# scale image based on rp radius criteria (for GZ)
        self.bg_image	    = single_image()		# add backgrounds (only possible for 5 SDSS bands at the moment)
#============ SET ORIGINAL IMAGE ======================#
        all_images,self.openlist  = sunpy.sunpy__load.load_all_broadband_images(filename,camera=camera,openlist=openlist,use_nonscatter=use_nonscatter)

        #to_nu                     = ((self.lambda_eff**2 ) / (speedoflight_m)) #* pixel_area_in_str
        to_nu = (self.ewidth_lambda/self.ewidth_nu)
        to_microjanskies          = (1.0e6) * to_nu * (1.0e26)                 # 1 muJy/str (1Jy = 1e-26 W/m^2/Hz)

        this_image = all_images[band,:,:]
        this_image = this_image * to_microjanskies 		# to microjanskies / str

        if True:	#verbose:
            print "SUNRISE calculated the abmag for this system to be: {:.2f}".format(self.filter_data.AB_mag_nonscatter0[band])

        self.sunrise_image.init_image(this_image, self, comoving_to_phys_fov=False)
            # assume now that all images are in micro-Janskies per str

        self.convolve_with_psf(add_psf=add_psf)

        #self.add_gaussian_psf(add_psf=add_psf)  add_gaussian_psf now called in convolve_with_psf, if appropriate
        self.rebin_to_physical_scale(rebin_phys=rebin_phys)
        self.add_noise(add_noise=add_noise, sn_limit=sn_limit, sky_sig=sky_sig)
        self.calc_r_petro(r_petro_kpc=r_petro_kpc, resize_rp=resize_rp)
        self.resize_image_from_rp(resize_rp=resize_rp)

        self.seed = seed
        self.bg_failed= False
        self.seed = self.add_background(seed=self.seed, add_background=add_background, rebin_gz=rebin_gz, n_target_pixels=n_target_pixels, fix_seed=fix_seed)

        end_time   = time.time()
        #print "init images + adding realism took "+str(end_time - start_time)+" seconds"
        num_label = len(bb_label)

        if verbose:
            print "preparing to save "+filename[:filename.index(bb_label)]+'synthetic_image_'+filename[filename.index(bb_label)+num_label:filename.index('.fits')]+'_band_'+str(self.band)+'_camera_'+str(camera)+'_'+str(int(self.seed))+'.fits'

        if save_fits:
            if custom_fitsfile != None:
                self.save_bgimage_fits_mujyas(custom_fitsfile,add_noise=add_noise,add_background=add_background)
                self.fitsfn = custom_fitsfile
            else:
                orig_dir=filename[:filename.index('broadband')]
                if bg_tag!=None:
                    outputfitsfile = orig_dir+output_label+'synthetic_image_'+filename[filename.index(bb_label)+num_label:filename.index('.fits')]+'_band_'+str(self.band)+'_camera_'+str(camera)+'_bg_'+str(int(bg_tag))+'.fits'
                else:
                    outputfitsfile = orig_dir+output_label+'synthetic_image_'+filename[filename.index(bb_label)+num_label:filename.index('.fits')]+'_band_'+str(self.band)+'_camera_'+str(camera)+'_bg_'+str(int(self.seed))+'.fits'
                self.save_bgimage_fits(outputfitsfile)
                self.fitsfn=outputfitsfile
            del self.sunrise_image, self.psf_image, self.rebinned_image, self.noisy_image, self.nmag_image, self.rp_image
            gc.collect()


    def convolve_with_psf(self, add_psf=True):
        if add_psf:
            if self.telescope.psf_fits_file != None:
                #first, rebin to psf pixel scale

                n_pixel_orig = self.sunrise_image.n_pixels
                n_pixel_new = self.sunrise_image.n_pixels*self.sunrise_image.pixel_in_arcsec/self.telescope.psf_pixsize_arcsec

                #print np.sum(self.sunrise_image.image)
                new_image = congrid(self.sunrise_image.image, (n_pixel_new,n_pixel_new))
                #print np.sum(new_image)

                #second, convolve with PSF
                if CONVOLVE_TYPE=='astropy':
                    #astropy.convolution.convolve()
                    print "convolving with astropy"
                    conv_im = convolve_fft(new_image,self.telescope.psf_kernel,boundary='fill',fill_value=0.0,normalize_kernel=True)  #boundary option?
                    #print np.sum(conv_im)
                else:
                    #scipy.signal.convolve2d()
                    conv_im = convolve(new_image,self.telescope.psf_kernel/np.sum(self.telescope.psf_kernel),boundary='fill',fillvalue=0.0,mode='same')  #boundary option?

                self.psf_image.init_image(conv_im,self)
                del new_image, conv_im
            else:
                self.add_gaussian_psf(add_psf=add_psf)
        else:
            self.psf_image.init_image(self.sunrise_image.image, self)
        gc.collect()


    def add_gaussian_psf(self, add_psf=True, sample_factor=1.0):		# operates on sunrise_image -> creates psf_image
        if add_psf:
            current_psf_sigma_pixels = self.telescope.psf_fwhm_arcsec * (1.0/2.355) / self.sunrise_image.pixel_in_arcsec

            if current_psf_sigma_pixels<8:	# want the psf sigma to be resolved with (at least) 8 pixels...
                target_psf_sigma_pixels  = 8.0
                n_pixel_new = np.floor(self.sunrise_image.n_pixels * target_psf_sigma_pixels / current_psf_sigma_pixels )

                if n_pixel_new > 1500:		# for speed, beyond this, the PSF is already very small...
                    n_pixel_new = 1500

                target_psf_sigma_pixels = n_pixel_new * current_psf_sigma_pixels / self.sunrise_image.n_pixels

                new_image = congrid(self.sunrise_image.image,  (n_pixel_new, n_pixel_new) )
                current_psf_sigma_pixels = target_psf_sigma_pixels * (
                            (self.sunrise_image.n_pixels * target_psf_sigma_pixels
                             / current_psf_sigma_pixels) / n_pixel_new )
            else:
                    new_image = self.sunrise_image.image

            psf_image = np.zeros_like( new_image ) * 1.0
            dummy = sp.ndimage.filters.gaussian_filter(new_image,
            current_psf_sigma_pixels, output=psf_image, mode='constant')

            self.psf_image.init_image(psf_image, self)
            del new_image, psf_image, dummy
        else:
            self.psf_image.init_image(self.sunrise_image.image, self)
        gc.collect()


    def rebin_to_physical_scale(self, rebin_phys=True):
	if rebin_phys:
	    n_pixel_new = np.floor( ( self.psf_image.pixel_in_arcsec / self.telescope.pixelsize_arcsec )  * self.psf_image.n_pixels )
	    rebinned_image = congrid(self.psf_image.image,  (n_pixel_new, n_pixel_new) )
  	    self.rebinned_image.init_image(rebinned_image, self) 
            del n_pixel_new, rebinned_image
            gc.collect()
	else:
	    self.rebinned_image.init_image(self.psf_image.image, self)
       

    def add_noise(self, add_noise=True, sky_sig=None, sn_limit=25.0):
	if add_noise:
	    if sky_sig==None:
	        total_flux 	= np.sum( self.rebinned_image.image )
	        area 		= 1.0 * self.rebinned_image.n_pixels * self.rebinned_image.n_pixels
	        sky_sig 	= np.sqrt( (total_flux / sn_limit)**2 / (area**2 ) )

	    noise_image 	=  sky_sig * np.random.randn( self.rebinned_image.n_pixels, self.rebinned_image.n_pixels ) 
	    new_image = self.rebinned_image.image + noise_image
	    self.noisy_image.init_image(new_image, self)
            del noise_image, new_image
            gc.collect()
	else:
	    self.noisy_image.init_image(self.rebinned_image.image, self)




    def calc_r_petro(self, r_petro_kpc=None, resize_rp=True):		# rename to "set_r_petro"
        " this routine is not working well.  Must manually set r_p until this is fixed..."
        if ( resize_rp==False ):
            r_petro_kpc = 1.0
        elif(r_petro_kpc==None):
            #RadiusObject 	= RadialInfo(self.noisy_image.n_pixels, self.noisy_image.image)
            r_petro_kpc = RadiusObject.PetroRadius * self.noisy_image.pixel_in_kpc    # do this outside of the RadialInfo class'
            if verbose:
                print " we've calculated a r_p of "+str(r_petro_kpc)
            del RadiusObject
            gc.collect()
        else:
            r_petro_kpc = r_petro_kpc
            if r_petro_kpc < 3.0:
                r_petro_kpc = 3.0
            if r_petro_kpc > 100.0:
                r_petro_kpc = 100.0


        r_petro_pixels = r_petro_kpc / self.noisy_image.pixel_in_kpc

        self.r_petro_pixels = r_petro_pixels
        self.r_petro_kpc    = r_petro_kpc


    def resize_image_from_rp(self, resize_rp=True, resize_factor=0.016, max_rp=100.0):
        if resize_rp:
            if self.r_petro_kpc < max_rp:
                rp_pixel_in_kpc = resize_factor * self.r_petro_kpc	# The target scale; was 0.008, upping to 0.016 for GZ based on feedback
            else:
                self.r_petro_kpc = max_rp
                rp_pixel_in_kpc = resize_factor * self.r_petro_kpc
            
            Ntotal_new = int( (self.noisy_image.pixel_in_kpc / rp_pixel_in_kpc ) * self.noisy_image.n_pixels )
            rebinned_image = congrid(self.noisy_image.image            ,  (Ntotal_new, Ntotal_new) )

            diff = n_pixels_galaxy_zoo - Ntotal_new		#
            if diff >= 0:
                shift = int(np.floor(1.0*diff/2.0))
        
                lp = shift
                up = shift + Ntotal_new
                tmp_image = np.zeros( (n_pixels_galaxy_zoo, n_pixels_galaxy_zoo) )
                tmp_image[lp:up,lp:up] = rebinned_image[0:Ntotal_new, 0:Ntotal_new]
                rp_image = tmp_image
            else:
                shift = int( np.floor(-1.0*diff/2.0) )
                lp = int(shift)
                up = int(shift+n_pixels_galaxy_zoo)
                rp_image = rebinned_image[lp:up, lp:up]


            self.rp_image.init_image(rp_image, self, fov = (1.0*n_pixels_galaxy_zoo)*(resize_factor * self.r_petro_kpc) )
            del rebinned_image, rp_image
            gc.collect()
        else:
            self.rp_image.init_image(self.noisy_image.image, self, fov=self.noisy_image.pixel_in_kpc*self.noisy_image.n_pixels)

	
    def add_background(self, seed=1, add_background=True, rebin_gz=False, n_target_pixels=n_pixels_galaxy_zoo, fix_seed=True):
        if add_background and (len(backgrounds[self.band]) > 0):
            bg_image = 10.0*self.rp_image.image		# dummy values for while loop condition
        
            tot_bg = np.sum(bg_image)
            tot_img= np.sum(self.rp_image.image)
            tol_fac = 1.0
        
            while(tot_bg > tol_fac*tot_img):
            #=== load *full* bg image, and its properties ===#
                bg_filename = (backgrounds[self.band])[0]
                if not (os.path.isfile(bg_filename)):
                    print "  Background files were not found...  "
                    print "  The standard files used in Torrey al. (2015), Snyder et al., (2015) and Genel et al., (2014) ..."
                    print "  can be downloaded using the download_backgrounds routine or manually from:  "
                    print "     http://illustris.rc.fas.harvard.edu/data/illustris_images_aux/backgrounds/SDSS_backgrounds/J113959.99+300000.0-u.fits "
                    print "     http://illustris.rc.fas.harvard.edu/data/illustris_images_aux/backgrounds/SDSS_backgrounds/J113959.99+300000.0-g.fits "
                    print "  "
            
                file = fits.open(bg_filename) ;       # was pyfits.open(bg_filename) ;
                header = file[0].header ;
                pixsize = get_pixelsize_arcsec(header) ;
                Nx = header.get('NAXIS2') ; Ny = header.get('NAXIS1')
                    
                    #=== figure out how much of the image to extract ===#
                Npix_get = np.floor(self.rp_image.n_pixels * self.rp_image.pixel_in_arcsec / pixsize)
  
                im = file[0].data 	# this is in some native units (nmaggies, for SDSS )
                halfval_i = np.floor(np.float(Nx)/1.3)
                halfval_j = np.floor(np.float(Ny)/1.3)
                np.random.seed(seed=int(seed))
                starti = np.random.random_integers(5,halfval_i)
                startj = np.random.random_integers(5,halfval_j)
                bg_image_raw = im[starti:starti+Npix_get,startj:startj+Npix_get]        # the extracted patch...
                #=== need to convert to microJy / str ===#
                bg_image_muJy = bg_image_raw * 10.0**(-0.4*(bg_zpt[self.band_name][0]- 23.9 ))       # if you got your zero points right, this is now in muJy
                pixel_area_in_str       = pixsize**2 / n_arcsec_per_str
                bg_image = bg_image_muJy / pixel_area_in_str
                #=== need to rebin bg_image  ===#
                bg_image = congrid(bg_image, (self.rp_image.n_pixels, self.rp_image.n_pixels))

                #=== compare sum(bg_image) to sum(self.rp_image.image) ===#
                if (fix_seed):
                    tot_bg = 0      # if seed is fixed, no need for brightness check...
                else:
                    tot_bg = np.sum(bg_image)
                    tot_img= np.sum(self.rp_image.image)
                    if(tot_bg > tol_fac*tot_img):
                        seed+=1
		
            new_image = bg_image + self.rp_image.image
            new_image[ new_image < self.rp_image.image.min() ] = self.rp_image.image.min()
            if (new_image.mean() > (5*self.rp_image.image.mean()) ):
                self.bg_failed=True
            del im, bg_image_raw, bg_image_muJy

        else:
            new_image = self.rp_image.image


        if rebin_gz:
            new_image = congrid( new_image, (n_target_pixels, n_target_pixels) )

        self.bg_image.init_image(new_image, self, fov = self.rp_image.pixel_in_kpc * self.rp_image.n_pixels)
        del new_image
        gc.collect()
        return seed



    def save_bgimage_fits(self,outputfitsfile, save_img_in_muJy=False):
	""" Written by G. Snyder 8/4/2014 to output FITS files from Sunpy module """
        theobj = self.bg_image
        image = np.copy( theobj.return_image() )		# in muJy / str

        pixel_area_in_str = theobj.pixel_in_arcsec**2 / n_arcsec_per_str
        image *= pixel_area_in_str      # in muJy
	print np.sum(image)
        if save_img_in_muJy == False:
	    print bg_zpt[self.band_name]
            if len(bg_zpt[self.band_name]) > 0:
                image = image / ( 10.0**(-0.4*(bg_zpt[self.band_name][0]- 23.9 )) ) 
		print ( 10.0**(-0.4*(bg_zpt[self.band_name][0]- 23.9 )) )
        else:
            print 'saving image in muJy!!!!!'
        print " "
        print " "
	print image.shape
        print np.sum(image) 
#        print 22.5 - 2.5*np.log10( np.sum(image) )
#        print -2.5*np.log10( np.sum(image) )
        print " "

        primhdu = fits.PrimaryHDU(image) ; primhdu.header.update('IMUNIT','NMAGGIE',comment='approx 3.63e-6 Jy')
        primhdu.header.update('ABABSZP',22.5,'For Final Image')  #THIS SHOULD BE CORRECT FOR NANOMAGGIE IMAGES ONLY
#       primhdu.header.update('ORIGZP',theobj.ab_abs_zeropoint,'For Original Image')
        primhdu.header.update('PIXSCALE',theobj.pixel_in_arcsec,'For Final Image, arcsec')
        primhdu.header.update('PIXORIG', theobj.camera_pixel_in_arcsec, 'For Original Image, arcsec')
        primhdu.header.update('PIXKPC',theobj.pixel_in_kpc, 'KPC')
        primhdu.header.update('ORIGKPC',self.sunrise_image.pixel_in_kpc,'For Original Image, KPC')
        primhdu.header.update('NPIX',theobj.n_pixels)
        primhdu.header.update('NPIXORIG',self.sunrise_image.n_pixels)

        primhdu.header.update('REDSHIFT',self.cosmology.redshift)
        primhdu.header.update('LUMDIST' ,self.cosmology.lum_dist, 'MPC')
        primhdu.header.update('ANGDIST' ,self.cosmology.ang_diam_dist, 'MPC')
        primhdu.header.update('PSCALE'  ,self.cosmology.kpc_per_arcsec,'KPC')
    
        primhdu.header.update('H0',self.cosmology.H0)
        primhdu.header.update('WM',self.cosmology.WM)
        primhdu.header.update('WV',self.cosmology.WV)

        primhdu.header.update('PSFFWHM',self.telescope.psf_fwhm_arcsec,'arcsec')
        primhdu.header.update('TPIX',self.telescope.pixelsize_arcsec,'arcsec')

        primhdu.header.update('FILTER', self.band_name)
        primhdu.header.update('FILE',self.filename)
#        primhdu.update_ext_name('SYNTHETIC_IMAGE')
	primhdu.header.name="SYNTHETIC_IMAGE"

#	primhdu.header[keyword] = value

        #Optionally, we can save additional images alongside these final ones
        #e.g., the raw sunrise image below
        #simhdu = pyfits.ImageHDU(self.sunriseimage, header=self.image_header) ; zhdu.update_ext_name('SIMULATED_IMAGE')
        #newlist = pyfits.HDUList([primhdu, simhdu])

        #create HDU List container

        newlist = fits.HDUList([primhdu])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            #save container to file, overwriting as needed
            newlist.writeto(outputfitsfile,clobber=True)


    def save_bgimage_fits_mujyas(self,outputfitsfile, save_img_in_muJy=False,add_noise=False, add_background=False):
	""" Written by G. Snyder 8/4/2014 to output FITS files from Sunpy module """
        """ Updated 9/24/2015 """

        theobj = self.bg_image
        image = np.copy( theobj.return_image() )		# in muJy / str
        image *= 1.0/n_arcsec_per_str      # in muJy/Arcsec**2

        sunobj = self.sunrise_image
        sunimage = np.copy(sunobj.return_image() )
        sunimage *= 1.0/n_arcsec_per_str

        #print theobj.pixel_in_arcsec
        AB_zeropoint = -2.5*np.log10(theobj.pixel_in_arcsec**2) - 2.5*(-6.0) + 2.5*np.log10(3631.0)  #for image in muJy/Arcsec**2
        total_apparent_mag = -2.5*np.log10(np.sum(image)) + AB_zeropoint
        total_absolute_mag = -2.5*np.log10(np.sum(image)) + AB_zeropoint - self.cosmology.distance_modulus
        sunrise_absolute_mag = self.sunrise_absolute_mag

        sun_AB_app_zp = -2.5*np.log10(sunobj.pixel_in_arcsec**2) - 2.5*(-6.0) + 2.5*np.log10(3631.0)
        sun_AB_cam_zp = -2.5*np.log10(sunobj.camera_pixel_in_arcsec**2) - 2.5*(-6.0) + 2.5*np.log10(3631.0)
        sun_AB_abs_zp = sun_AB_app_zp - self.cosmology.distance_modulus

        sunrise_image_camera_mag = -2.5*np.log10(np.sum(sunimage)) + sun_AB_cam_zp
        sunrise_image_apparent_mag = -2.5*np.log10(np.sum(sunimage)) + sun_AB_app_zp
        sunrise_image_absolute_mag = -2.5*np.log10(np.sum(sunimage)) + sun_AB_abs_zp

        primhdu = fits.PrimaryHDU(np.float32(image)) ; primhdu.header.update('IMUNIT','muJy/SqArcsec',comment='microjanskies per square arcsecond')
        primhdu.header.update('ABZP',round(AB_zeropoint,6),'For Final Image')
        primhdu.header.update('PIXSCALE',round(theobj.pixel_in_arcsec,6),'For Final Image, arcsec')
        primhdu.header.update('PIXORIG', round(theobj.camera_pixel_in_arcsec,6), 'For Original Image, arcsec')
        primhdu.header.update('PIXKPC',round(theobj.pixel_in_kpc,6), 'KPC')
        primhdu.header.update('ORIGKPC',round(self.sunrise_image.pixel_in_kpc,6),'For Original Image, KPC')
        primhdu.header.update('NPIX',theobj.n_pixels)
        primhdu.header.update('NPIXORIG',self.sunrise_image.n_pixels)

        primhdu.header.update('REDSHIFT',self.cosmology.redshift)
        primhdu.header.update('LUMDIST' ,round(self.cosmology.lum_dist,6), 'MPC')
        primhdu.header.update('ANGDIST' ,round(self.cosmology.ang_diam_dist,6), 'MPC')
        primhdu.header.update('PSCALE'  ,round(self.cosmology.kpc_per_arcsec,6),'KPC')
        primhdu.header.update('DISTMOD'  ,round(self.cosmology.distance_modulus,6),'Mag')

        primhdu.header.update('H0',round(self.cosmology.H0,6))
        primhdu.header.update('WM',round(self.cosmology.WM,6))
        primhdu.header.update('WV',round(self.cosmology.WV,6))

        if self.telescope.psf_fits_file==None:
            primhdu.header.update('PSFFWHM',round(self.telescope.psf_fwhm_arcsec,6),'arcsec')
        else:
            primhdu.header.update('PSFFILE',os.path.join(os.path.basename(os.path.dirname(self.telescope.psf_fits_file)),os.path.basename(self.telescope.psf_fits_file)))

        primhdu.header.update('TPIX',round(self.telescope.pixelsize_arcsec,6),'arcsec')

        primhdu.header.update('FILTER', self.band_name)
        primhdu.header.update('FILE',self.filename)
        primhdu.header.update('EFLAMBDA',round(self.lambda_eff*1.0e6,6),'filter effective wavelength [microns]')

        primhdu.header.update('MAG', round(total_apparent_mag,6), 'AB system')
        primhdu.header.update('ABSMAG', round(total_absolute_mag,6), 'AB system')
        primhdu.header.update('SUNMAG', round(sunrise_absolute_mag,6), 'from spectrum, Note: excludes Lyman absorption')
        primhdu.header.update('SUNCMAG', round(sunrise_image_camera_mag,6), 'from image, camera mag')
        primhdu.header.update('SUNAPMAG', round(sunrise_image_apparent_mag,6), 'from image, apparent mag')
        primhdu.header.update('SUABSMAG', round(sunrise_image_absolute_mag,6), 'from image, absolute mag')

        if add_noise==False and add_background==False:
            primhdu.header.update('SKYSIG', 0.0, 'image units')
        elif sky_sig != None:
            primhdu.header.update('SKYSIG', round(self.sky_sig,6), 'image units')

        if add_background==True:
            primhdu.header.update('BGFILE', os.path.basename(backgrounds[self.band]))

        camera_param_cards = self.param_header.cards[13:]
        for card in camera_param_cards:
            #print card
            primhdu.header.append(card)


        primhdu.update_ext_name('SYNTHETIC_IMAGE')

        #Optionally, we can save additional images alongside these final ones
        #e.g., the raw sunrise image below
        #simhdu = pyfits.ImageHDU(self.sunriseimage, header=self.image_header) ; zhdu.update_ext_name('SIMULATED_IMAGE')
        #newlist = pyfits.HDUList([primhdu, simhdu])

        #create HDU List container

        newlist = fits.HDUList([primhdu])

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            #save container to file, overwriting as needed
            newlist.writeto(outputfitsfile,clobber=True)


def get_pixelsize_arcsec(header):
    cd1_1 = header.get('CD1_1')  # come in degrees	
    cd1_2 = header.get('CD1_2')

    if cd1_2==None:
        cd1_2 = header.get('CD2_2')

    try:
        pix_arcsec = 3600.0*(cd1_1**2 + cd1_2**2)**0.5
    except:
        print "WARNING!!! SETTING PIXEL SCALE MANUALLY!"
        pix_arcsec = 0.05

    return pix_arcsec



from scipy.optimize import curve_fit

def my_fit(r, a, b, c):
     return a * np.exp(-r / b) + c







class RadialInfo:
    """ Class for giving radial profile info for rp calcultions """
    def __init__(self,N, image, num_pts=100, max_pixels_for_fit=100000):
        
        self.Npix = N
        self.RadiusGrid = np.linspace(0.01,1.0*N,num=num_pts)
        self.PetroRatio = np.ones_like( self.RadiusGrid)

        xgrid = np.linspace(float(-self.Npix)/2.0 + 0.5,float(self.Npix)/2.0 - 0.5,num=self.Npix)
        xsquare = np.zeros((self.Npix,self.Npix))
        ysquare = np.zeros_like(xsquare)
        ones = np.ones((self.Npix,self.Npix))
        for j in range(self.Npix):
            xsquare[j,:] = xgrid
            ysquare[:,j] = xgrid
    
        self.rsquare = (xsquare**2 + ysquare**2)**0.5

        x0 = np.array(self.rsquare).flatten()
        y0 = np.array(image).flatten()
        
        #print x0.shape
        x0 = x0[ y0 > 0 ]
        y0 = y0[ y0 > 0 ]
        
        if x0.shape[0] > max_pixels_for_fit:
            index_list = np.arange( x0.shape[0] )
            index_list = np.random.choice( index_list, max_pixels_for_fit)
        
            x0 = x0[index_list]
            y0 = y0[index_list]
        
        popt, pcov = curve_fit(my_fit, x0, np.log10(y0))
        y1 = 10.0**(my_fit(self.RadiusGrid, *popt))
        fake_image1 = 10.0**(my_fit(self.rsquare, *popt))
        y2 = 10.0**(my_fit(self.RadiusGrid, *popt)) - 10.0**popt[2]
        fake_image2 = 10.0**(my_fit(self.rsquare, *popt)) - 10.0**popt[2]
        
        
        y1sum = np.zeros_like(y1)
        y2sum = np.zeros_like(y2)
        for index,val in enumerate(y1[:-1]):
            if index==0:
                this_r = 0.5*(self.RadiusGrid[index]+self.RadiusGrid[index+1])
                y1sum[index] = 3.14159 * this_r**2 * y1[index]
                y2sum[index] = 3.14159 * this_r**2 * y2[index]

            else:
                y1sum[index] = y1sum[index-1] + 3.14159 * (self.RadiusGrid[index+1]**2 - self.RadiusGrid[index]**2 ) * (y1[index] + y1[index+1])/2.0
                y2sum[index] = y2sum[index-1] + 3.14159 * (self.RadiusGrid[index+1]**2 - self.RadiusGrid[index]**2 ) * (y2[index] + y2[index+1])/2.0


        for index,val in enumerate(y1sum[:-1]):
            this_r = (0.5*(self.RadiusGrid[index]+self.RadiusGrid[index+1]))
            y1sum[index] = val / (3.14159 * this_r **2 )
            y2sum[index] = y2sum[index] / (3.14159 * this_r **2 )
        
        
        self.PetroRatio = np.array(y2/y2sum)
        self.PetroRatio[np.isnan(self.PetroRatio)] = 0.0
        self.PetroRatio[np.isinf(self.PetroRatio)] = 0.0
        self.Pind        = np.argmin( np.absolute( np.flipud(self.PetroRatio) - 0.2) )
        self.PetroRadius = np.flipud(self.RadiusGrid)[self.Pind]
        
        if verbose:
            print y2
            print " Saving Figure ..."
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(22,5))
            ax = fig.add_subplot(1,5,1)
            print x0.shape
            ax.plot(x0, y0, 'ro', ms=5)
            ax.plot(self.RadiusGrid, y1, 'g', lw=3, ls='-')
            ax.plot(self.RadiusGrid, y2, 'b', lw=3, ls='-')
            ax.plot(self.RadiusGrid, y1sum, 'g', lw=3, ls='-.')
            ax.plot(self.RadiusGrid, y2sum, 'b', lw=3, ls='-.')
            ax.plot([self.PetroRadius, self.PetroRadius], [1,1e20], 'k', lw=1, ls='-')
            ax.set_yscale('log')
            ax.set_ylim([1e7,1e13])
            ax = fig.add_subplot(1,5,2)
            ax.plot(self.RadiusGrid, y1/y1sum, 'g', lw=3)
            ax.plot(self.RadiusGrid, y2/y2sum, 'b', lw=3)
            ax.plot([self.PetroRadius, self.PetroRadius], [-10,10], 'k', lw=1, ls='-')
            ax.plot(self.RadiusGrid, np.ones_like(self.RadiusGrid) * 0.2 )
            ax.set_ylim([0,1])
            ax = fig.add_subplot(1,5,3)
            ax.imshow( np.log10(image), vmin=7, vmax=13 )
            ax = fig.add_subplot(1,5,4)
            ax.imshow( np.log10(fake_image1), vmin=7, vmax=13  )
            ax = fig.add_subplot(1,5,5)
            ax.imshow( np.log10(fake_image2), vmin=7, vmax=13  )
            fig.savefig('temp1.png')
            fig.clf()
            plt.close()
#            print " Figure has been saved ... "




        
#        for i,rad in enumerate(self.RadiusGrid):
#            tf_annulus = np.logical_and( self.rsquare < 1.25*rad, self.rsquare > 0.80*rad )
#            tf_annulus = np.logical_and( tf_annulus, image > min_img_thresh )
#            self.annulus_indices.append(   np.where(tf_annulus) ) #np.logical_and( self.rsquare < 1.25*rad, self.rsquare > 0.80*rad )) )#

#            tf_int = np.logical_and( self.rsquare < rad, image > min_img_thresh )
#            self.interior_indices.append(  np.where(tf_int) )
#            self.annulus_sums.append(  np.sum(ones[self.annulus_indices[i]] ) )
#            self.interior_sums.append( np.sum(ones[self.interior_indices[i]]) )
#            this_sum = np.sum( image[self.interior_indices[i]])



#        for radius in RadiusObject.RadiusGrid:
#            pflux_annulus 	= image[ self.annulus_indices[i]  ]
#            pflux_interior 	= image[ self.interior_indices[i] ]
#            self.sumI_r[i] 	= np.sum(pflux_interior)
#            if(self.annulus_sums[i]*self.interior_sums[i] != 0.0):
#                self.AnnulusSB[i]  = (np.sum(pflux_annulus)/self.annulus_sums[i])
#                self.IntSB[i]      = (np.sum(pflux_interior)/self.interior_sums[i])
#                self.PetroRatio[i] = (np.sum(pflux_annulus)/self.annulus_sums[i])/(np.sum(pflux_interior)/self.interior_sums[i])








class fits_header:
  def __init__(self, filename):
    if (not os.path.exists(filename)):
      print "file not found:", filename
      sys.exit()

    hdulist = fits.open(filename)
    self.info = hdulist.info()



def my_fits_open(filename):
    if (not os.path.exists(filename)):
        print "file not found:", filename
        sys.exit()

    return fits.open(filename)










#============ COSMOLOGY PARAMETERS =====================#
# cosmology class:
#
#  used to track (i) the cosmological parameters and 
#  (ii) image properties set by our adopted cosmology
#
#  This class is used to distinguish features of the telescope
#  (e.g., pixel size in arcseconds) from features of our 
# adopted cosmology (e.g.,image kpc per arcsec)
#
#=======================================================#
class cosmology:
    def __init__(self, redshift, H0=70.4, WM=0.2726, WV=0.7274):
        self.H0=H0
        self.WM=WM
        self.WV=WV
        self.redshift = redshift
        self.lum_dist       = (cosmocalc.cosmocalc(self.redshift, H0=self.H0, WM=self.WM, WV=self.WV))['DL_Mpc']          ## luminosity dist in mpc
        self.ang_diam_dist  = (cosmocalc.cosmocalc(self.redshift, H0=self.H0, WM=self.WM, WV=self.WV))['DA_Mpc']          ## 
        self.kpc_per_arcsec = (cosmocalc.cosmocalc(self.redshift, H0=self.H0, WM=self.WM, WV=self.WV))['PS_kpc']
        self.distance_modulus = 5.0 * ( np.log10(self.lum_dist*1.0e6) - 1.0 )



#============ TELESCOPE PARAMETERS =====================#
# telescope class:
#
# used to track the psf size in arcsec and pixelsize in arcsec
#=======================================================#
class telescope:
    def __init__(self, psf_fwhm_arcsec, pixelsize_arcsec, psf_fits, psf_pixsize_arcsec,rebin_phys,add_psf,psf_truncate_pixels,psf_hdu_num):
        self.psf_fwhm_arcsec  = psf_fwhm_arcsec
        self.pixelsize_arcsec = pixelsize_arcsec
        self.psf_truncate_pixels = psf_truncate_pixels
        self.psf_hdu_num = psf_hdu_num

        self.psf_fits_file = None
        self.psf_kernel = None
        self.psf_pixsize_arcsec = None

        #future upgrade:  pass kernel directly and pixel scale instead?
        if psf_fits != None:
            self.psf_fits_file = psf_fits
            orig_psf_hdu = fits.open(psf_fits,ignore_missing_end=True)[psf_hdu_num]                
            orig_psf_kernel = orig_psf_hdu.data
            #some psfs come in cubes... what does this parameter mean?  for STDPSF: fiducial detector positions... want near center
            if orig_psf_kernel.ndim==3:
                npsf = orig_psf_kernel.shape[0]
                orig_psf_kernel = orig_psf_kernel[npsf/2,:,:]

            if psf_truncate_pixels != None:
                psfc = orig_psf_kernel.shape[0]/2
                st = self.psf_truncate_pixels
                orig_psf_kernel = orig_psf_kernel[psfc-st:psfc+st,psfc-st:psfc+st]

            #psf kernel shape must be odd for astropy.convolve??
            if orig_psf_kernel.shape[0] % 2 == 0:
                new_psf_shape = orig_psf_kernel.shape[0]-1
                self.psf_kernel = congrid(orig_psf_kernel,(new_psf_shape,new_psf_shape))
            else:
                self.psf_kernel = orig_psf_kernel

            assert( self.psf_kernel.shape[0] % 2 != 0)
            assert (psf_pixsize_arcsec != None)
            self.psf_pixsize_arcsec = psf_pixsize_arcsec

            if (self.psf_pixsize_arcsec > self.pixelsize_arcsec) and (rebin_phys==True) and (add_psf==True):
                print "WARNING: you are requesting to rebin an image to a higher resolution than the requested PSF file supports. OK if this is desired behavior."



#=====================================================#
# single_image class:
# 
# This class is used to host and track the properties for 
# a single image (one galaxy, one band, one level of realism).
# This class tracks important image traits, such as the 
# image array itself, the field of view, number of pixels,
# ab_zeropoint, pixel scale, etc.
# 
# When new images are created (e.g., when psf bluring is 
# done on the original image) a new "single_image" instance 
# is created.  
#
# The synthetic_image class (defined below) contains 
# several instances of this single_image class
#
#=====================================================#
class single_image:
    def __init__(self):
        self.image_exists = False

    def init_image(self, image, parent_obj, fov=None, comoving_to_phys_fov=False):
        self.image              = image
        self.n_pixels           = image.shape[0]
        if fov==None:
	    if comoving_to_phys_fov:
                self.pixel_in_kpc           = parent_obj.param_header.get('linear_fov') / self.n_pixels / (parent_obj.cosmology.redshift+1)
	    else:
		self.pixel_in_kpc           = parent_obj.param_header.get('linear_fov') / self.n_pixels
        else:
            self.pixel_in_kpc           = fov / self.n_pixels
        self.pixel_in_arcsec    = self.pixel_in_kpc / parent_obj.cosmology.kpc_per_arcsec
        self.image_exists       = True
        self.camera_pixel_in_arcsec = (self.pixel_in_kpc / parent_obj.param_header.get('cameradist') ) * 2.06e5

	pixel_in_sr = (1e3*self.pixel_in_kpc /10.0)**2
	image_in_muJy =  self.image  * pixel_in_sr		# should now have muJy
        tot_img_in_Jy = np.sum(image_in_muJy) / 1e6		# now have total image flux in Jy
	abmag = -2.5 * np.log10(tot_img_in_Jy / 3631 )
#	print "the ab magnitude of this image is :"+str(abmag)


    def calc_ab_abs_zero(self, parent_obj):
        lambda_eff_in_m         = parent_obj.lambda_eff
        pixel_area_in_str       = self.camera_pixel_in_arcsec**2 / n_arcsec_per_str
        cameradist_in_kpc       = parent_obj.param_header.get('cameradist')

        to_nu                     = ((lambda_eff_in_m**2 ) / (speedoflight_m))* pixel_area_in_str
        to_microjanskies          = (1.0e6) * to_nu * (1.0e26)             # 1 Jy = 1e-26 W/m^2/Hz
        to_microjanskies_at_10pc  = to_microjanskies * (cameradist_in_kpc / abs_dist)**2  

        ab_abs_zeropoint = 23.90 - (2.5*np.log10(to_microjanskies_at_10pc))             
        self.ab_abs_zeropoint = ab_abs_zeropoint

    def convert_orig_to_nanomaggies(self, parent_obj):
        distance_factor = (10.0 / (parent_obj.cosmology.lum_dist * 1.0e6))**2
        orig_to_nmaggies = distance_factor * 10.0**(0.4*(22.5 - self.ab_abs_zeropoint) )
        self.image_in_nmaggies = self.image * orig_to_nmaggies

    def return_image(self):
        return self.image


def return_img_nanomaggies_to_orig(image_nm, lum_dist, ab_abs_zeropoint):
    distance_factor = (10.0 / (lum_dist * 1.0e6))**2
    orig_to_nmaggies = distance_factor * 10.0**(0.4*(22.5 - ab_abs_zeropoint) )
    return image_nm / orig_to_nmaggies






def congrid(a, newdims, centre=False, minusone=False):
    ''' Slimmed down version of congrid as originally obtained from:
		http://wiki.scipy.org/Cookbook/Rebinning
    '''
    if not a.dtype in [np.float64, np.float32]:
        a = np.cast[float](a)

    m1 = np.cast[int](minusone)
    ofs = np.cast[int](centre) * 0.5
    old = np.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
        return None
    newdims = np.asarray( newdims, dtype=float )
    dimlist = []



    for i in range( ndims ):
        base = np.arange( newdims[i] )
        dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
    # specify old dims
    olddims = [np.arange(i, dtype = np.float) for i in list( a.shape )]

    # first interpolation - for ndims = any
    mint = scipy.interpolate.interp1d( olddims[-1], a, kind='linear', bounds_error=False, fill_value=0.0 )
    newa = mint( dimlist[-1] )

    trorder = [ndims - 1] + range( ndims - 1 )
    for i in range( ndims - 2, -1, -1 ):
        newa = newa.transpose( trorder )

        mint = scipy.interpolate.interp1d( olddims[i], newa, kind='linear', bounds_error=False, fill_value=0.0 )
        newa = mint( dimlist[i] )

    if ndims > 1:
        # need one more transpose to return to original dimensions
        newa = newa.transpose( trorder )

    return newa

def download_backgrounds():
    if not os.path.exists('./data'):
        os.makedirs('./data')
    if not os.path.exists('./data/SDSS_backgrounds'):
        os.makedirs('./data/SDSS_backgrounds')
    if not os.path.exists('./data/HST_backgrounds'):
        os.makedirs('./data/HST_backgrounds')
    for this_background in backgrounds:
        if len(this_background) > 0:
            if not (os.path.isfile(this_background[0])):
                url=dl_base+this_background[0][len(bg_base):]
                this_file = wget.download(url)
                os.rename(this_file, this_background[0])

                print url
