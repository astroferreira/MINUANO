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

from astropy.convolution import convolve
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
import math

verbose = 0

def get_footprints(hdu):
    """Calculates positions of the corners of the science extensions of some image 'im_name' in sky space"""
    footprints = []
    

    hdr = hdu.header
    wcs = WCS(hdr)
    footprint = wcs.calc_footprint(hdr)
    footprints.append(footprint)

    return footprints

def bounds(footprint_list):
    """Calculate RA/Dec bounding box properties from multiple RA/Dec points"""
    
    # flatten list of extensions into numpy array of all corner positions
    merged = [ext for image in footprint_list for ext in image]
    merged = np.vstack(merged)
    ras, decs = merged.T
    
    # Compute width/height
    delta_ra = (max(ras)-min(ras))
    delta_dec = max(decs)-min(decs)
    # Compute midpoints
    ra_midpt = (max(ras)+min(ras))/2.
    dec_midpt = (max(decs)+min(decs))/2.
    
    return ra_midpt, dec_midpt, delta_ra, delta_dec

def get_error_mask(catalog, max_error):
    """Returns a mask for rows in catalog where RA and Dec error are less than max_error"""
    ra_mask = catalog['ra_error']< max_error
    dec_mask = catalog['dec_error'] < max_error
    #mag_mask = catalog['phot_g_mean_mag'] > min_mag
    mask = ra_mask & dec_mask
    #     print('Cutting sources with error higher than {}'.format(max_error))
    #     print('Number of sources befor filtering: {}\nAfter filtering: {}\n'.format(len(mask),sum(mask)))
    return mask


if __name__ == '__main__':

    SCIi = 1
    path_to_file = sys.argv[1]
    ref_cat = sys.argv[2]

    input_image = fits.open(path_to_file)
    input_data = input_image[SCIi].data.copy()
    wcs = WCS(input_image[1].header)
    input_wcs = HSTWCS(input_image, SCIi)


    print('Fix Stars statured pixels')
    from photutils.detection import DAOStarFinder
    from astropy.stats import sigma_clipped_stats
    from photutils.datasets import load_star_image

    mean, median, std = sigma_clipped_stats(input_data, sigma=3.0)  
    daofind = DAOStarFinder(fwhm=6.0, threshold=10. * std, roundlo=-0.05, roundhi=0.05)
    sources = daofind(input_data-median)
    
    RA, DEC = input_wcs.pixel_to_world_values(sources['xcentroid'], sources['ycentroid'])
    coords = SkyCoord(ra=RA, dec=DEC, unit=u.deg)
    cat = Table([RA, DEC, sources['xcentroid'], sources['ycentroid']])
    
    #cat = cat.to_table()
    #cat['RA']  = cat['sky_centroid'].ra.degree  
    #cat['DEC'] = cat['sky_centroid'].dec.degree 
    cat.rename_column('col0', 'RA')
    cat.rename_column('col1', 'DEC')
    cat.rename_column('xcentroid', 'x')
    cat.rename_column('ycentroid', 'y')
    print(cat)
    print(f'Found {len(cat)} sources')

    if ref_cat is None:
        print('Querying GAIA...')
        Gaia.ROW_LIMIT = 200
        Gaia.MAIN_GAIA_TABLE = "gaiaedr3.gaia_source"
        footprint_list = [get_footprints(input_image[SCIi])]
        ra_midpt, dec_midpt, delta_ra, delta_dec = bounds(footprint_list)
        coord = SkyCoord(ra=ra_midpt, dec=dec_midpt, unit=u.deg)
        
        width = Quantity(delta_ra, u.deg)
        height = Quantity(delta_dec, u.deg)
        r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
        ras = r['ra']
        decs = r['dec']

        mask = get_error_mask(r, 10)

        gaia_table = Table([ras[mask], decs[mask]], names=['RA', 'DEC']) 
        print(gaia_table)
        gaia_coords = SkyCoord(ra=ras, dec=decs, unit=u.deg)
    else:

        gaia_table = Table.read(ref_cat)
        print(gaia_table)
        gaia_coords =  SkyCoord(ra=gaia_table['RA'], dec=gaia_table['DEC'], unit=u.deg)

    
    print(f'Found {len(gaia_table)} GAIA sources...')

    #match catalogs

    match = TPMatch(searchrad=3, separation=0.06, tolerance=10, use2dhist=False, xoffset=0.05, yoffset=0.05)
    input_wcs_corrector = FITSWCS(input_wcs)
    ridx, iidx = match(gaia_table, cat, input_wcs_corrector)
    print('Number of matches:', len(ridx), len(iidx))

    seps = []
    for ri, ii in zip(ridx, iidx):
        sep = gaia_coords[ri].separation(coords[ii])
        seps.append(sep.to(u.arcsec).value)
    
    
    seps = np.array(seps)
    print(seps, seps.mean())
    mask_sep = seps < 30
    print('Applying rejection of separations > 10 pixels')
    print(f'Reducing selection to {len(cat[iidx][mask_sep])}')
    print(f'Mean OFFSET: {seps[mask_sep].mean()}')
    
    
    fitted = fit_wcs(gaia_table[ridx][mask_sep], cat[iidx][mask_sep], input_wcs_corrector, sigma=0.06, fitgeom='general')
    aligned_imwcs = fitted.wcs
    
    print(fitted)
    print(fitted.meta['shift'])
    print(fitted.meta['matrix'])
    if 1:
        if verbose:
            print('OLD WCS')
            print(input_wcs)
            print('----------')
            print('NEW WCS')
            print(aligned_imwcs)
            
        updatehdr.update_wcs(input_image, SCIi, aligned_imwcs, wcsname='TWEAK', reusename=True, verbose=True)
        input_image.writeto(path_to_file, overwrite=True)
        print(f'Finished aligning {path_to_file} to GAIA.')
