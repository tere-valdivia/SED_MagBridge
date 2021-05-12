from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy import units as u
from astropy.nddata import Cutout2D

folder = '/Users/terevaldivia/Documents/uchile/Cursos/Investigacion Dirigida/Maryland/Maryland/HERITAGE/'

filenames = ['SMC.HERITAGE.PACS100.img', 'SMC.HERITAGE.PACS160.img', 'SMC_250um_combined_20121215_img',
             'SMC_350um_combined_20121215_img', 'SMC_500um_combined_20121215_img']
# ra_E = 31.7229859 * u.deg
# dec_E = -74.7577409 * u.deg
ra_F = 33.6594937 * u.deg
dec_F = -74.3604975 * u.deg
# MagBridgeE_pos = SkyCoord(ra_E, dec_E, frame='fk5')
MagBridgeF_pos = SkyCoord(ra_F, dec_F, frame='fk5')
beamdict_key = ['100', '160', '250', '350', '500']
beamdict_Herschel = {'100': ((9*u.arcsec).to(u.deg)).value, '160': ((14*u.arcsec).to(u.deg)).value, '250': (
    (22*u.arcsec).to(u.deg)).value, '350': ((30*u.arcsec).to(u.deg)).value, '500': ((43*u.arcsec).to(u.deg)).value}

for filename, key in zip(filenames, beamdict_key):
    data = fits.getdata(folder+filename+'.fits')
    header = fits.getheader(folder+filename+'.fits')
    radius = 6 * u.arcmin
    wcs = WCS(header)

    cutout = Cutout2D(data, MagBridgeF_pos, (radius*2, radius*2), wcs=wcs)
    cutout = Cutout2D(data, MagBridgeF_pos, (radius*2, radius*2), wcs=wcs)
    newdata = cutout.data
    newheader = cutout.wcs.to_header()
    # newheader['OBJECT'] = 'MagBri_E'
    newheader['OBJECT'] = 'MagBri_F'
    newheader['BMAJ'] = beamdict_Herschel[key]
    newheader['BMIN'] = beamdict_Herschel[key]
    newheader['BPA'] = 0.0

    #
    # if 'OBJECT' in header:
    #     newheader['OBJECT'] = header['OBJECT']
    if 'BUNIT' in header:
        newheader['BUNIT'] = header['BUNIT']
    elif 'ZUNITS' in header:
        newheader['ZUNITS'] = header['ZUNITS']
    if 'JANSCALE' in header:
        newheader['JANSCALE'] = header['JANSCALE']
    if 'WAVELNTH' in header:
        newheader['WAVELNTH'] = header['WAVELNTH']
    # cutout_filename = 'images/MagBridgeE_' + filename + '.fits'
    cutout_filename = 'images/MagBridgeF_' + filename + '.fits'
    hdunew = fits.PrimaryHDU(data=newdata, header=newheader)
    hdunew.writeto(cutout_filename, overwrite=True)
