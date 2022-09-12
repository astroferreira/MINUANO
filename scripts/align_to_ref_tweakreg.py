__version__ = 0.1
__author__  = 'Nathan Adams & Leonardo Ferreira'
__email__   = 'leonardo.ferreira@nottingham.ac.uk'
__date__    = '2022 07 16'

"""
    Align JWST WCS to GAIA

    Solution based on Dan Coe's notebook: source here
    For debugging use the notebook
"""

import os
import sys
import numpy as np

import tweakwcs
from tweakwcs import fit_wcs, FITSWCS, TPMatch
from stwcs.wcsutil import HSTWCS

from drizzlepac import updatehdr

from astropy.convolution import convolve_fft as convolve
from astropy.convolution import Gaussian2DKernel

from photutils.segmentation import detect_sources
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import SourceCatalog

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.table import Table
from astropy.units import Quantity
from drizzlepac import tweakreg
import math


if __name__ == '__main__':

    path_to_file = sys.argv[1]
    path_to_file_ref = sys.argv[2]
    
    input_image = fits.open(path_to_file)
    input_data = input_image[1].data
    input_wcs = HSTWCS(input_image, 1)

    ref_image = fits.open(path_to_file_ref)
    ref_data = ref_image[0].data
    print(ref_data)
    ref_wcs = HSTWCS(ref_image, 0)


#    kernel = Gaussian2DKernel(x_stddev=3)
#    print('Fix Stars statured pixels input')
#    input_data[input_data==0] = np.nan
#    input_data_fixed = convolve(input_data, kernel)

#    print('Fix Stars statured pixels ref')
#    ref_data[ref_data==0] = np.nan
#    ref_data_fixed = convolve(ref_data, kernel)
    

    bkg_estimator = MedianBackground()
    bkg = Background2D(ref_data, (50, 50), filter_size=(3, 3),
                    bkg_estimator=bkg_estimator)

    threshold = 5 * bkg.background_rms

    print('Detecting Souces to Align...')
    segment_map = detect_sources(ref_data, threshold, npixels=20)
    ref_cat = SourceCatalog(ref_data, segment_map, wcs=ref_wcs)
    ref_cat = ref_cat.to_table()
    ref_cat['RA']  = ref_cat['sky_centroid'].ra.degree  
    ref_cat['DEC'] = ref_cat['sky_centroid'].dec.degree 
    ref_cat.rename_column('xcentroid', 'x')
    ref_cat.rename_column('ycentroid', 'y')

    tbl = Table([ref_cat['RA'], ref_cat['DEC']])
    tbl.write('tweakreg_ref.cat', format='ascii.fast_commented_header', overwrite=True)
    print(f'Found {len(ref_cat)} sources')

    print(tbl)

    tweakreg.TweakReg(path_to_file, # Pass input images
                  updatehdr=True, # update header with new WCS solution
                  imagefindcfg={'threshold':10.,'conv_width':12.0},# Detection parameters, threshold varies for different data
                  separation=0., # Allow for very small shifts
                  searchrad = 1.0,
                  searchunits='arcseconds',
                  nclip=3,
                  sigma=3,
                  tolerance=10,
                  refcat='tweakreg_ref.cat', # Use user supplied catalog (Gaia)
                  clean=True, # Get rid of intermediate files
                  interactive=False,
                  residplot='both',
                  see2dplot=True,
                  shiftfile=True, # Save out shift file (so we can look at shifts later)
                  wcsname='RELICS', # Give our WCS a new name
                  reusename=True,
                  verbose=True,
                  fitgeometry='general') # Use the 6 parameter fit


    input_image = fits.open(path_to_file)
    from reproject import reproject_interp
    array, footprint = reproject_interp(input_image[1], ref_image[0].header)
#    pixel_aligned_file = input_image.copy()

#   pixel_aligned_file[1].data = array

    
    fits.writeto(f'{path_to_file.replace("1overf", "1overf_sci")}', array, ref_image[0].header, overwrite=True)

    print(f'Finished aligning {path_to_file} to {path_to_file_ref}.')
