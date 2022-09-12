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
from reproject import reproject_interp

import astropy.units as u
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astropy.table import Table
from astropy.units import Quantity
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
    bkg = Background2D(input_data, (50, 50), filter_size=(3, 3),
                    bkg_estimator=bkg_estimator)

    threshold = 10 * bkg.background_rms

    print('Detecting Souces to Align...')
    segment_map = detect_sources(input_data, threshold, npixels=50)
    input_cat = SourceCatalog(input_data, segment_map, wcs=input_wcs)
    input_cat = input_cat.to_table()
    input_cat['RA']  = input_cat['sky_centroid'].ra.degree  
    input_cat['DEC'] = input_cat['sky_centroid'].dec.degree 
    input_cat.rename_column('xcentroid', 'x')
    input_cat.rename_column('ycentroid', 'y')

    print(f'Found {len(input_cat)} sources')

    bkg_estimator = MedianBackground()
    bkg = Background2D(ref_data, (50, 50), filter_size=(3, 3),
                    bkg_estimator=bkg_estimator)

    threshold = 10 * bkg.background_rms

    print('Detecting Souces to Align...')
    segment_map = detect_sources(ref_data, threshold, npixels=30)
    ref_cat = SourceCatalog(ref_data, segment_map, wcs=ref_wcs)
    ref_cat = ref_cat.to_table()
    ref_cat['RA']  = ref_cat['sky_centroid'].ra.degree  
    ref_cat['DEC'] = ref_cat['sky_centroid'].dec.degree 
    ref_cat.rename_column('xcentroid', 'x')
    ref_cat.rename_column('ycentroid', 'y')

    print(f'Found {len(ref_cat)} sources')

    #match catalogs
    match = TPMatch(searchrad=1, separation=0.05, tolerance=5, use2dhist=False, xoffset=0.1, yoffset=0.1)
    input_wcs_corrector = FITSWCS(input_wcs)
    ridx, iidx = match(ref_cat, input_cat, input_wcs_corrector)
    print('Number of matches:', len(ridx), len(iidx))

    seps = []
    for ri, ii in zip(ridx, iidx):
        sep = ref_cat[ri]['sky_centroid'].separation(input_cat[ii]['sky_centroid'])
        seps.append(sep.to(u.arcsec).value/0.03)
    
    
    seps = np.array(seps)
    print(seps, seps.mean())
    mask_sep = seps < 100
    print('Applying rejection of separations > 10 pixels')
    print(f'Reducing selection to {len(input_cat[iidx][mask_sep])}')
    print(f'Mean OFFSET: {seps[mask_sep].mean()}')
    
    
    if 1:
        aligned_imwcs = fit_wcs(ref_cat[ridx][mask_sep], input_cat[iidx][mask_sep], input_wcs_corrector).wcs
    
        print('OLD WCS')
        print(input_wcs)
        print('----------')
        print('NEW WCS')
        print(aligned_imwcs)
        
        updatehdr.update_wcs(input_image, 1, aligned_imwcs, wcsname='TWEAK', reusename=True, verbose=True)
        input_image.close()


        print(f'Finished aligning {path_to_file} to {path_to_file_ref}.')
