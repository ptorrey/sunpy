#!/usr/bin/env python
"""  A very simple examle script to download, open, and plot an image
     from the Illustris Image host site.

     Author:  P. Torrey (ptorrey@mit.edu) 1/24/15
"""

import numpy as np
import os
import wget					# https://pypi.python.org/pypi/wget -- can be installed via PIP
import sunpy.sunpy__load as sunpy__load		# 
import sunpy.sunpy__plot as sunpy__plot



dl_base='http://illustris.rc.fas.harvard.edu/data/'		# the base data directory

try:
    catalog = np.loadtxt('directory_catalog_135.txt',
		    dtype={'names'  : ('subdirs', 'galaxy_numbers', 'galaxy_masses'),
                           'formats': ('S3', 'i10', 'f8')})
except:
    wget.download(dl_base+"illustris_images_aux/directory_catalog_135.txt")
    catalog = np.loadtxt('directory_catalog_135.txt',
		    dtype={'names'  : ('subdirs', 'galaxy_numbers', 'galaxy_masses'),
	                   'formats': ('S3', 'i10','f8')})


all_subdirs = catalog['subdirs']
all_galnrs  = catalog['galaxy_numbers']

for index,galnr in enumerate(all_galnrs[:1]):
    url=dl_base+"illustris_images/subdir_"+str(all_subdirs[index])+"/broadband_"+str(galnr)+".fits"

    print "Want to load galaxy ID="+str(galnr)+" from folder "+str(all_subdirs[index])
    print "  path="+url
    print " "
    
    if( !(os.path.isfile("./broadband_"+str(galnr)+".fits")) )
        filename = wget.download(url)
    else:
        filename = "./broadband_"+str(galnr)+".fits"

    # retrieve the image.  Could be used for non-parametric fitting, plotting, etc.
    image = sunpy__load.load_broadband_image(filename,band='B_Johnson.res')

    # plot an SDSS gri band idealized image
    sunpy__plot.plot_sdss_gri(filename, savefile='./sdss_gri_'+str(galnr)+'.png')

    os.remove("./"+filename)

