#sunpy

Module for opening, manipulating, and plotting mock galaxy images produced with SUNRISE.
This module has been developed to support the Illustris Image Pipeline (http://arxiv.org/abs/1411.3717).  


## Table of Contents 
 - [Data Access](#data-access)
 - [Example Usage](#example-usage)
 - [Contributors](#contributors)
 

## Data Access
The sunpy module operates on FITS files produced by the SUNRISE radiative transfer code.
The module was developed with the specific intent of being used in the Illustris Image Pipeline 
publically available data set.  
Access to the data can be obtained through the Illustris Website (http://www.illustris-project.org/) 
by registering under the "Galaxy Observatory" tab.  Once this is done, users are able to search and 
download individual and batches of FITS files.

Alternatively, the mock images can also be downloaded directly from URL's that take the following form
```
http://illustris.rc.fas.harvard.edu/data/illustris_images/subdir_543/broadband_12345.fits
```
In this path there are two relevant numbers:  the subdir number (always 3 digits, in this case 543) and the galaxy number (in this case, 12345).  
In order for you to download a FITS file directly, you'll need to know both numbers.  
An ASCII file with three columns: subdir number, galaxy number, and the log of the stellar mass can be found here
```
http://illustris.rc.fas.harvard.edu/data/illustris_images_aux/directory_catalog_135.txt
```
Once you have the list of corresponding subdir numbers and galaxy numbers, it is straight forward to piece together the full URL from which the content can be downloaded.  
This can executed in practice either through your browser, or using wget or a similar command.


## Example Usage
Example scripts showing how to access, open, manipulate, and plot the data can be found in this repository under the examples folder.



## Contributors

Contributors:
- Paul Torrey (ptorrey@mit.edu)
- Greg Snyder (gsnyder@stsci.edu)



Greg Snyder is responsible for the "image realism" components of this module including:
  - PSF blurring
  - Rebinning images to new pixel scale given input redshift and telescope angular resolution
  - Adding of random noise
  - Adding of background images (SDSS)


Paul Torrey is responsible for most of the FITS file loading, manipulation, and plotting routines.






