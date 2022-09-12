from jwst.datamodels import ImageModel, FlatModel, dqflags

from glob import glob
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import simple_norm,ImageNormalize
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint
from astropy.stats import SigmaClip
from photutils.background import Background2D, SExtractorBackground, MedianBackground

from scipy.ndimage import rotate


import os
import sys


def collapse_image(im, mask, dimension='y', sig=2.):
    """collapse an image along one dimension to check for striping.
    By default, collapse columns to show horizontal striping (collapsing
    along columns). Switch to vertical striping (collapsing along rows)
    with dimension='x' 
    Striping is measured as a sigma-clipped median of all unmasked pixels 
    in the row or column.
    Args:
        im (float array): image data array
        mask (bool array): image mask array, True where pixels should be 
            masked from the fit (where DQ>0, source flux has been masked, etc.)
        dimension (Optional [str]): specifies which dimension along which 
            to collapse the image. If 'y', collapses along columns to 
            measure horizontal striping. If 'x', collapses along rows to 
            measure vertical striping. Default is 'y'
        sig (Optional [float]): sigma to use in sigma clipping
    """
    # axis=1 results in array along y
    # axis=0 results in array along x

    if dimension == 'y':
#        collapsed = np.median(im, axis=1)
        res = sigma_clipped_stats(im, mask=mask, sigma=sig, cenfunc='median',
                                        stdfunc='std', axis=1)

    elif dimension == 'x':
#        collapsed = np.median(im, axis=0)
        res = sigma_clipped_stats(im, mask=mask, sigma=sig, cenfunc='median',
                                        stdfunc='std', axis=0)

    return res[1]
    








filename = sys.argv[1]
boxsize = int(sys.argv[2])
align = sys.argv[3]




image = fits.open(filename)


ang = 0
if align == '--align':
    if image['PRIMARY'].header['MODULE'] == 'A':
        ang = 30
    elif image['PRIMARY'].header['MODULE'] == 'B':
        ang = 29.9        

    data = rotate(image['SCI'].data, ang)
    DQ_rot = rotate(image['DQ'].data, ang)
    mask = np.zeros(data.shape, dtype=bool)
    mask[DQ_rot > 0] = True
    coveragemask = data == 0
    mask[coveragemask] = True
    slowaxis = abs(image['PRIMARY'].header['SLOWAXIS'])

    stripes=np.zeros(data.shape)
    divisor = 1
    loamp = np.linspace(0, data.shape[0]-data.shape[0] / divisor, divisor).astype(int)
    hiamp = np.linspace(data.shape[0] / divisor, data.shape[0], divisor).astype(int)
    numamps=loamp.size

    sigma_clip = SigmaClip(sigma=3, maxiters=5)
    threshold = detect_threshold(data, nsigma=10.0, sigma_clip=sigma_clip)
    segment_img = detect_sources(data, threshold, npixels=10)
    footprint = circular_footprint(radius=25)
    sourcemask = segment_img.make_source_mask(footprint=footprint)
    mask[sourcemask] = True


    for k in range(numamps):
        if slowaxis==1:
            maskedbksubdata = np.ma.array(data[loamp[k]:hiamp[k],:], mask=mask[loamp[k]:hiamp[k],:])
            stripes[loamp[k]:hiamp[k],:] = np.ma.median(maskedbksubdata, axis=(slowaxis-1),keepdims=True)
        elif slowaxis==2:
            maskedbksubdata = np.ma.array(data[:,loamp[k]:hiamp[k]], mask=mask[:,loamp[k]:hiamp[k]])
            stripes[:,loamp[k]:hiamp[k]] = np.ma.median(maskedbksubdata, axis=(slowaxis-1),keepdims=True)


    #data[np.isnan(data)] = 0
    stripes = rotate(stripes, -ang)[1902-1024:1902+1024, 1902-1024:1902+1024]

    data = image['SCI'].data.copy()
    data -= stripes
    sigma_clip = SigmaClip(sigma=3, maxiters=5)
    threshold = detect_threshold(data, nsigma=10.0, sigma_clip=sigma_clip)
    segment_img = detect_sources(data, threshold, npixels=10)
    footprint = circular_footprint(radius=25)
    sourcemask = segment_img.make_source_mask(footprint=footprint)
    bkg_estimator = SExtractorBackground()
    bkg = Background2D(data, (boxsize, boxsize), mask=sourcemask,  coverage_mask=image[1].data == 0, filter_size=(3, 3),
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)

    image[1].data = data-bkg.background
    image.writeto(filename.replace('_bg', '_bg_flat'), overwrite=True)
else:

    print(filename)
    model = ImageModel(filename)
    mask = np.zeros(model.data.shape, dtype=bool)
    mask[model.dq > 0] = True
    coveragemask = model.data == 0
    mask[coveragemask] = True

    #find sources
    sigma_clip = SigmaClip(sigma=3, maxiters=5)
    threshold = detect_threshold(model.data, nsigma=10.0, sigma_clip=sigma_clip)
    segment_img = detect_sources(model.data, threshold, npixels=10)
    footprint = circular_footprint(radius=25)
    sourcemask = segment_img.make_source_mask(footprint=footprint)

    data = model.data.copy()
    #remove striping
    horizontal_striping = collapse_image(data, mask, dimension='y')
    # remove horizontal striping, requires taking transpose of image
    temp_image = data.T.copy()
    temp_image -= horizontal_striping
    # transpose back
    data = temp_image.T

    # fit vertical striping, collapsing along rows
    vertical_striping = collapse_image(data, mask, dimension='x')
    data -= vertical_striping
    data[coveragemask] = 0


    bkg_estimator = SExtractorBackground()
    bkg = Background2D(data, (boxsize, boxsize), mask=sourcemask,  coverage_mask=coveragemask, filter_size=(3, 3),
                    sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)


    image[1].data = data-bkg.background
    image.writeto(filename.replace('_bg', '_bg_flat'), overwrite=True)