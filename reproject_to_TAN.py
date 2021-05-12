import numpy as np
from reproject import reproject_interp
from astropy.io import fits

hdu = fits.open('images/MagBri_F_sm18.fits')[0]
new_header = hdu.header.copy()
new_header['CTYPE1'] = 'RA---TAN'
new_header['CTYPE2'] = 'DEC--TAN'

new_image, footprint = reproject_interp(hdu, new_header)

fits.writeto('images/MagBri_F_sm18_reprojected.fits',
             new_image, new_header)
