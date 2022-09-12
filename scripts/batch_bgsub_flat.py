
import numpy as np
from astropy.io import fits
from matplotlib import pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import simple_norm,ImageNormalize
from astropy.stats import sigma_clipped_stats, SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint


import sys

filename = sys.argv[1]
image = fits.open(filename)


print(f'running for {filename}')

sigma_clip = SigmaClip(sigma=3.0, maxiters=10)
threshold = detect_threshold(image[1].data, nsigma=2.0, sigma_clip=sigma_clip)
segment_img = detect_sources(image[1].data, threshold, npixels=10)
footprint = circular_footprint(radius=10)
mask = segment_img.make_source_mask(footprint=footprint)
mean, median, std = sigma_clipped_stats(image[1].data, sigma=3.0, mask=mask)
image[1].data -= median
image[1].header['BGSUB'] = median
print(f'Subtracting BG[{median}] from {filename}')
image.writeto(filename.replace('_cal_1overf', '_cal_1overf_bg'), overwrite=True)