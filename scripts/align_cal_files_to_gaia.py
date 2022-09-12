from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
from astropy.units import Quantity
from matplotlib import pyplot as plt
import numpy as np
from photutils.detection import DAOStarFinder
from astropy.stats import sigma_clipped_stats
from photutils.datasets import load_star_image

from astropy.visualization.mpl_normalize import simple_norm,ImageNormalize
from astropy.visualization import SqrtStretch
from astropy.visualization import simple_norm
from astroquery.gaia import Gaia

from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel

from photutils.segmentation import detect_sources
from photutils.background import Background2D, MedianBackground
from photutils.segmentation import SourceCatalog

import tweakwcs
from tweakwcs import fit_wcs, FITSWCS, TPMatch
from stwcs.wcsutil import HSTWCS
from drizzlepac import updatehdr

from astroquery.gaia import Gaia


from photutils.segmentation import SourceFinder


def query_gaia(input_image, SCIi=0):
    wcs = WCS(input_image[SCIi].header)
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

    mask = get_error_mask(r, 30)

    gaia_table = Table([ras[mask], decs[mask]], names=['RA', 'DEC']) 
    print(gaia_table)
    gaia_coords = SkyCoord(ra=ras[mask], dec=decs[mask], unit=u.deg)
    return gaia_table, gaia_coords

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
    #mag_mask = catalog['phot_g_mean_mag'] > 16
    mask = ra_mask & dec_mask# & mag_mask
    #     print('Cutting sources with error higher than {}'.format(max_error))
    #     print('Number of sources befor filtering: {}\nAfter filtering: {}\n'.format(len(mask),sum(mask)))
    return mask

def reconstruct_saturated_stars(hdu, use_max=True):

    kernel = Gaussian2DKernel(x_stddev=3)

    
    data = hdu['SCI'].data
    
    if use_max:
        DQ = hdu['DQ'].data
        max_data = data[~np.isnan(data)].max()
        data[DQ == 2] = max_data/10
    else:
        data[data==0] = np.nan
        

    return convolve(data, kernel)


def detect_round_sources(data):

    mean, median, std = sigma_clipped_stats(data, sigma=3.0)  
    daofind = DAOStarFinder(fwhm=10.0, threshold=10. * std, roundlo=-0.5, roundhi=0.5)
    sources = daofind(data-median)

    return sources

from photutils.segmentation import SourceFinder
from photutils.centroids import centroid_sources, centroid_com
def detect_all_sources(data, dq=None):
    
    kernel = Gaussian2DKernel(x_stddev=1)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                    bkg_estimator=bkg_estimator)

    threshold = 20 * bkg.background_rms

    print('Detecting Souces to Align...')
    
    if dq is not None:
        data[dq>0] = np.nan
        data[data<=0] = np.nan
    
    #data = convolve(data, kernel)
     
    finder = SourceFinder(npixels=100, progress_bar=False)
    segment_map = finder(data-bkg.background, threshold)
    #segment_map = detect_sources(data, threshold, npixels=50)
    sources = SourceCatalog(data, segment_map)
    sources = sources.to_table()

    sources = sources[np.where(~np.isnan(sources['xcentroid'].value))]
    
    x, y = centroid_sources(data, sources['xcentroid'], sources['ycentroid'], boxsize=41, centroid_func=centroid_com)
    
    sources['xcentroid'] = x
    sources['ycentroid'] = y
    
    return sources


def detect_all_sources2(data, dq=None):
        
    
    
    kernel = Gaussian2DKernel(x_stddev=1)
    bkg_estimator = MedianBackground()
    bkg = Background2D(data, (50, 50), filter_size=(3, 3),
                    bkg_estimator=bkg_estimator)

    threshold = 20 * bkg.background_rms

    print('Detecting Souces to Align...')
    
    if dq is not None:
        data[dq>0] = np.nan
        data[data<=0] = np.nan
    
    #data = convolve(data, kernel)
     
    finder = SourceFinder(npixels=100, progress_bar=False)
    segment_map = finder(data-bkg.background, threshold)
    #segment_map = detect_sources(data, threshold, npixels=50)
    sources = SourceCatalog(data, segment_map)
    sources = sources.to_table()
    
    return sources

import sys

if __name__ == '__main__':

    filename = sys.argv[1]
    sep_max = float(sys.argv[2])
    tolerance = int(sys.argv[3])
    refcat = None
    
    try:
        refcat = sys.argv[4]
    except:
        pass

    print(refcat)

    SCIi = 1
    pixel_scale = 0.03
    verbose=0
    hdu = fits.open(filename)
    
    wcs = HSTWCS(hdu, SCIi)
    #wcs = WCS(hdu[SCIi].header)

    data = hdu[1].data#reconstruct_saturated_stars(hdu, use_max=refcat is None)

    try:
        sources = detect_all_sources(data, dq=hdu['DQ'].data)#detect_round_sources(data)
    except:
        sources = detect_all_sources(data)#detect_round_sources(data)


    RA, DEC = wcs.pixel_to_world_values(sources['xcentroid'], sources['ycentroid'])
    coords = SkyCoord(ra=RA, dec=DEC, unit=u.deg)
    cat = Table([RA, DEC, sources['xcentroid'], sources['ycentroid']])
    cat.rename_column('col0', 'RA')
    cat.rename_column('col1', 'DEC')
    cat.rename_column('xcentroid', 'x')
    cat.rename_column('ycentroid', 'y')

    print(f'Found {len(cat)} sources in the image')

    if refcat is None:
        gaia_table, gaia_coords = query_gaia(hdu, SCIi=1)
        gaia_x, gaia_y = wcs.world_to_pixel(gaia_coords)
    else:
        gaia_table = Table.read(refcat)
        gaia_coords = SkyCoord(ra=gaia_table['RA'], dec=gaia_table['DEC'], unit=u.deg)
        #gaia_x = gaia_table['x']
        #gaia_y = gaia_table['y']


    #units in PIXELS
    match = TPMatch(searchrad=5, separation=10, tolerance=tolerance, use2dhist=False)
    input_wcs_corrector = FITSWCS(wcs)
    ridx, iidx = match(gaia_table, cat, input_wcs_corrector)

    print('Number of matches:', len(ridx), len(iidx))
   
    seps = []
    for ri, ii in zip(ridx, iidx):
        sep = gaia_coords[ri].separation(coords[ii])
        seps.append(sep.to(u.arcsec).value)

    seps = np.array(seps)
    mask = seps < sep_max
    if len(seps[seps<sep_max]) > 2:
        info = fit_wcs(gaia_table[ridx][mask], cat[iidx][mask],
         input_wcs_corrector, fitgeom='general', nclip=1)
        print(f'rmse: {info.meta["fit_info"]["rmse"]}')
        updatehdr.update_wcs(filename, 1, info.wcs, wcsname='TWEAK', reusename=True, verbose=True)
    else:
        with open(filename.replace('fits', 'shift'), 'w') as shift_file:
            shift_file.write('MEAN: rmse: inf\n')
        print('Not enough sources to use as alignment reference, check the shift file')
        exit()
    



    if verbose:
        f = plt.figure(figsize=(10, 10))
        norm = ImageNormalize(stretch=SqrtStretch(), vmin=-0.00, vmax=0.5)
        plt.imshow(data, norm=norm, cmap='gray_r')
        plt.scatter(cat['x'][iidx], cat['y'][iidx], marker='x' , color='r', alpha=0.7)
        plt.scatter(gaia_x[ridx], gaia_y[ridx], marker='s',  facecolors='none', color='cyan', alpha=0.7)
        plt.savefig('align_debug.png')


    
    with open(filename.replace('fits', 'shift'), 'w') as shift_file:
        for sep in seps[mask]:
            shift_file.write(f'{sep}\n')

        shift_file.write(f'[{len(seps[mask])}] MEAN:{np.median(seps) : .4} - {np.std(seps) : .4} - {hdu[0].header["FILTER"]}-{hdu[0].header["DETECTOR"]}\n')
            
        shift_file.write('SHIFT\n')
        shift_file.write(f'DX: {info.meta["shift"][0]}, DY {info.meta["shift"][1]}\n')
        shift_file.write(f'rmse: {info.meta["fit_info"]["rmse"]}')
        
        shift_file.flush()

