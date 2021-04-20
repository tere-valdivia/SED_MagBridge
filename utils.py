import numpy as np
import astropy.units as u
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D


def get_beam(header):
    bmaj = header['bmaj'] * u.deg
    bmin = header['bmin'] * u.deg
    bpa = header['bpa'] * u.deg
    beamarea = 1.133 * bmaj * bmin
    return bmaj, bmin, bpa, beamarea


def convert_to_Jy_beam(header, initunit):
    bmaj, bmin, bpa, beamarea = get_beam(header)
    if initunit == 'MJy/sr':
        return (1 * u.MJy/u.sr).to(u.Jy/u.beam, equivalencies=u.beam_angular_area(beamarea)).value
    if initunit == 'Jy/beam':
        return 1


def convert_to_Jy(header, initunit):
    bmaj, bmin, bpa, beamarea = get_beam(header)
    if initunit == 'MJy/sr':
        brightness_conversion = (1 * u.MJy/u.sr).to(u.Jy/u.beam,
                                                    equivalencies=u.beam_angular_area(beamarea)).value
        return brightness_conversion * ((header['CDELT2']*u.deg)**2/beamarea)
    if initunit == 'Jy/beam':
        return ((header['CDELT2']*u.deg)**2/beamarea)


def cut_image_with_beam(filename, cropcenter, cropwidth, cropheight, overwrite=False):
    # Open the data and header
    data = fits.getdata(filename+'.fits')
    header = fits.getheader(filename+'.fits')
    wcs = WCS(header)
    # define the cutout region
    ra, dec = cropcenter
    position = SkyCoord(ra, dec)
    width = int(abs(cropwidth.to(u.deg).value / header['CDELT1']))
    height = int(abs(cropheight.to(u.deg).value / header['CDELT2']))
    size = (width, height)
    # do the cutout
    cutout = Cutout2D(data, position, size, wcs=wcs)
    # create the new data and header
    newdata = cutout.data
    newheader = cutout.wcs.to_header()
    # fill the header with the beam and brightness units
    newheader['bmaj'] = header['bmaj']
    newheader['bmin'] = header['bmin']
    newheader['bpa'] = header['bpa']
    newheader['bunit'] = header['bunit']
    newheader['btype'] = header['btype']
    newheader['bscale'] = header['bscale']
    newheader['bzero'] = header['bzero']
    # save the new Data
    hdu = fits.PrimaryHDU(data=newdata, header=newheader)
    hdu.writeto(filename+'_cut_'+str(width)+'_'+str(height)+'.fits', overwrite=overwrite)


# filename = 'images/MagBridgeA_BoA_iter_RM_10_reprojected_43arcsec'
# center = np.array([25.96189167, -74.53989194]) * u.deg
# width = height = 6 * u.arcmin
# cut_image_with_beam(filename, center, width, height, overwrite=True)
