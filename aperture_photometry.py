import numpy as np
import pandas as pd
from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from utils import *
from ChauvenetRMS import *
from photutils import aperture_photometry, CircularAperture, SkyCircularAperture, SkyCircularAnnulus

"""
Be careful that the images have the beam in them
The rms will not be correctly estimated in the case there is an extended
primary beam correction in the map
Final flux and errors are in Jy

TODO: correct to take rms in ring
"""
####
# input parameters

folder = 'images/'
# imagelist = ['MagBridgeA_100um_43arcsec', 'MagBridgeA_160um_43arcsec', 'MagBridgeA_250um_43arcsec', 'MagBridgeA_350um_43arcsec',
#              'MagBridgeA_500um', 'MagBridgeA_BoA_iter_RM_10_reprojected_43arcsec_cut_79_79']  # must all be with the same resolution
imagelist = ['MagBridgeF_SMC.HERITAGE.PACS100.img_43', 'MagBridgeF_SMC.HERITAGE.PACS160.img_43', 'MagBridgeF_SMC_250um_combined_20121215_img_43', 'MagBridgeF_SMC_350um_combined_20121215_img_43',
             'MagBridgeF_SMC_500um_combined_20121215_img', 'MagBridgeF_MagBri_F_sm18_reprojected_43']  # must all be with the same resolution

# percentage of error in the flux calibration for each image
fluxcal_error = [0.1, 0.2, 0.08, 0.08, 0.08, 0.2]
# MagBridgeA
# aperture_center = np.array([25.96189167, -74.53989194]) * u.deg
# aperture_radius = 50 * u.arcsec  # arcsec
# MagBridgeF
aperture_center = np.array([33.6623458, -74.3574931]) * u.deg
aperture_radius = 46 * u.arcsec  # arcsec: covers full S/N=3 contour
rmslist = [2.0, 5.0, 0.5, 0.3, 0.2, 0.005]
# for different skies, do a list
sky_center = aperture_center
# MagBridgeA
# sky_in = np.ones(np.shape(imagelist)) * 2.5 * u.arcmin
# sky_out = np.ones(np.shape(imagelist)) * 3 * u.arcmin
# MagBridgeF
sky_in = np.ones(np.shape(imagelist)) * 1.5 * u.arcmin
sky_out = np.ones(np.shape(imagelist)) * 2.5 * u.arcmin
savefile = 'aphot_MagBridgeF'
####

assert len(imagelist) == len(fluxcal_error) and len(
    imagelist) == len(sky_in) and len(imagelist) == len(sky_out)

columns = ['aperture_sum', 'aperture_npix', 'aperture_flux', 'sky_sum',
           'sky_mean', 'sky_npix', 'sky_flux', 'rms', 'flux', 'flux_error_cal', 'flux_error_noise', 'flux_error']
results = pd.DataFrame(columns=columns)

for i in range(len(imagelist)):
    filename = folder + imagelist[i]
    fluxcal = fluxcal_error[i]
    r_in = sky_in[i]
    r_out = sky_out[i]

    # open the image
    header = fits.getheader(filename+'.fits')
    data = fits.getdata(filename+'.fits')
    wcs = WCS(header).celestial

    # determine the beam and conversion factor
    bmaj, bmin, bpa, beamarea = get_beam(header)
    if not 'BUNIT' in header.keys():
        units = header['ZUNITS']
    else:
        units = header['BUNIT']
        if units == '':
            units = header['ZUNITS']
    pixsize = header['CDELT2'] * header['PC2_2'] * u.deg
    conversion = convert_to_Jy(header, units)
    print(units)
    # calculate the rms
    rms, noisedata = calculatenoise(data)
    if i == len(imagelist) - 1:
        rms = rmslist[i]
    results.loc[filename, 'rms'] = rms

    # define the apertures
    position = SkyCoord(aperture_center[0], aperture_center[1])
    aperture = SkyCircularAperture(position, r=aperture_radius)
    position_sky = SkyCoord(aperture_center[0], aperture_center[1])
    aperture_sky = SkyCircularAnnulus(position_sky, r_in=r_in, r_out=r_out)
    pix_aperture = aperture.to_pixel(wcs)
    pix_sky = aperture_sky.to_pixel(wcs)
    results.loc[filename, 'aperture_npix'] = pix_aperture.area()
    results.loc[filename, 'sky_npix'] = pix_sky.area()

    # do the aperture photometry in source
    phot_table_aperture = aperture_photometry(data, pix_aperture)
    results.loc[filename, 'aperture_sum'] = phot_table_aperture['aperture_sum'][0]
    results.loc[filename, 'aperture_flux'] = results.loc[filename, 'aperture_sum'] * conversion

    # do the aperture photometry in background
    phot_table_sky = aperture_photometry(data, pix_sky)
    results.loc[filename, 'sky_sum'] = phot_table_sky['aperture_sum'][0]
    results.loc[filename, 'sky_mean'] = phot_table_sky['aperture_sum'][0] / \
        results.loc[filename, 'sky_npix']
    results.loc[filename, 'sky_flux'] = results.loc[filename, 'sky_mean'] * \
        results.loc[filename, 'aperture_npix'] * conversion

    # do the subtraction
    results.loc[filename, 'flux'] = results.loc[filename,
                                                'aperture_flux'] - results.loc[filename, 'sky_flux']
    # flux calibration error
    results.loc[filename, 'flux_error_cal'] = results.loc[filename, 'flux'] * fluxcal

    # noise error
    rms_jy = rms * convert_to_Jy_beam(header, units)
    results.loc[filename, 'flux_error_noise'] = rms_jy * \
        np.sqrt(results.loc[filename, 'aperture_npix'] * (pixsize**2/beamarea).value)
    results.loc[filename, 'flux_error'] = np.sqrt(
        results.loc[filename, 'flux_error_noise']**2 + results.loc[filename, 'flux_error_cal']**2)

    # save in table
results.to_csv(savefile+'.csv')
