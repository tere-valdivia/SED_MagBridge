import numpy as np
import os

filenames = ['SMC.HERITAGE.PACS100.img', 'SMC.HERITAGE.PACS160.img', 'SMC_250um_combined_20121215_img',
             'SMC_350um_combined_20121215_img']
for filename in filenames:
    name = 'images/MagBridgeE_' + filename
    importfits(name+'.fits', imagename=name + '.image')
    ia.open(name + '.image')
    im2 = ia.convolve2d(outfile=name + '_43.image',
                        major='43arcsec', minor='43arcsec', pa='0deg', targetres=True, overwrite=True)
    ia.close()
    im2.done()
    exportfits(imagename=name+'_43.image', fitsimage=name + '_43.fits')
    os.system('rm -r {}'.format(name + '.image'))
    os.system('rm -r {}'.format(name + '_43.image'))
