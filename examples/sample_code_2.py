#!/usr/bin/env python
"""  An examle script to download, open, and plot an image
     from the Illustris Image host site.  This routine also 
     adds images realism to create mock SDSS images by co-adding 
     real SDSS backgrounds.

     Author:  P. Torrey (ptorrey@mit.edu) 1/24/15
"""

import numpy as np
import sunpy.sunpy__load as sunpy__load
import sunpy.sunpy__plot as sunpy__plot
import os

my_api = "enter your own api key here !!!"
dl_base='http://www.illustris-project.org'

try:
    catalog = np.loadtxt('directory_catalog_135.txt',
		    dtype={'names'  : ('subdirs', 'galaxy_numbers', 'galaxy_masses'),
                           'formats': ('S3', 'i10', 'f8')})
except:
    os.system("wget "+dl_base+"/files/directory_catalog_135.txt")
    catalog = np.loadtxt('directory_catalog_135.txt',
		    dtype={'names'  : ('subdirs', 'galaxy_numbers', 'galaxy_masses'),
	                   'formats': ('S3', 'i10','f8')})

all_subdirs = catalog['subdirs']
all_galnrs  = catalog['galaxy_numbers']


common_args = { 
                'add_background':       False,          # default option = False, turn on one-by-one
                'add_noise':            False,
                'add_psf':              False,
                'rebin_phys':           False,
                'resize_rp':            False,
                'rebin_gz':             True,           # always true, all pngs 424x424
                'scale_min':            1e-10,          # was 1e-4
                'lupton_alpha':         2e-12,          # was 1e-4
                'lupton_Q':             10,             # was ~100
                'pixelsize_arcsec':     0.24,
                'psf_fwhm_arcsec':      1.0,
                'sn_limit':             25.0,           # super low noise value, dominated by background
                'sky_sig':              None,           #
                'redshift':             0.05,           # 
                'b_fac':                1.1, 
                'g_fac':                1.0, 
                'r_fac':                0.9,
                'camera':               1,
                'seed_boost':           1.0,
                'save_fits':            True
                }


def restore_common_args():
    common_args['rebin_phys']           = False
    common_args['add_psf']              = False
    common_args['add_background']       = False        
    common_args['add_noise']            = False


for index,galnr in enumerate(all_galnrs[:1]):
    cmd = 'wget --content-disposition --header="API-Key: '+my_api+'" "'+dl_base+ \
        '/api/Illustris-1/snapshots/135/subhalos/'+str(galnr)+  \
        '/stellar_mocks/broadband.fits"'

    if( not (os.path.isfile("./broadband_"+str(galnr)+".fits")) ):
	print "trying to download"
        os.system(cmd)
	print "did it download?"

    filename = "./broadband_"+str(galnr)+".fits"


    sunpy__plot.plot_sdss_gri(filename, savefile='./sdss_gri_'+str(galnr)+'.png')
    sunpy__plot.plot_synthetic_sdss_gri(filename, savefile='./synthetic_0_sdss_gri_'+str(galnr)+'.png' , **common_args)

    common_args['rebin_phys'] = True
    sunpy__plot.plot_synthetic_sdss_gri(filename, savefile='./synthetic_1_sdss_gri_'+str(galnr)+'.png' , **common_args)

    common_args['add_psf'] = True
    sunpy__plot.plot_synthetic_sdss_gri(filename, savefile='./synthetic_2_sdss_gri_'+str(galnr)+'.png' , **common_args)

    common_args['add_background'] = True
    common_args['add_noise'] = True
    sunpy__plot.plot_synthetic_sdss_gri(filename, savefile='./synthetic_3_sdss_gri_'+str(galnr)+'.png' , **common_args)


