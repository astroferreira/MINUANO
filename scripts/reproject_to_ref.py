from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from reproject import reproject_interp
import glob
import copy
from multiprocessing import Pool


verbose = 0
SCIi = 1
slices_to_project = [1, 2, 3, 4, 5, 6, 7]
path = '/home/ppxlf2/data/JWST/CEERS_REAL/science'
files = sorted(glob.glob('/home/ppxlf2/data/JWST/CEERS_REAL/science/*o001*i2d.fits'))

#print(files)
if verbose:
    #reference is F444W
    print(f'Reference file is: {files[-1]}')
    
names = [f.split('/')[-1].split('.fits')[0] for f in files]
refname = names[1]
if verbose:
    from matplotlib import pyplot as plt
    from astropy.visualization import make_lupton_rgb, SqrtStretch, LogStretch, hist, simple_norm
    norm = simple_norm(to_reproject.data, 'sqrt', min_percent=1, max_percent=99)

    f = plt.figure(figsize=(20, 30))

    ax1 = plt.subplot(1,2,1, projection=WCS(reference.header))
    ax1.imshow(reference.data, norm=norm, origin='lower')
    ax1.coords['ra'].set_axislabel('Right Ascension')
    ax1.coords['dec'].set_axislabel('Declination')
    ax1.set_title('F444W')

    ax2 = plt.subplot(1,2,2, projection=WCS(to_reproject.header))
    ax2.imshow(to_reproject.data, norm=norm, origin='lower')
    ax2.set_title('F090W')
    
    
reference = fits.open(files[1])
def align(name):
    print(name)

    if name == refname:
        print(name, 'exiting')
        return

    to_reproject = fits.open(f'{path}/{name}.fits')
    pixel_aligned_file = copy.deepcopy(to_reproject)
    
    for slice in slices_to_project:
        to_reproject[slice].header = to_reproject[SCIi].header
        array, footprint = reproject_interp(to_reproject[slice], reference[SCIi].header)
        
        (NAXIS1, NAXIS2) = array.shape
        print(NAXIS1, NAXIS2)
        pixel_aligned_file[slice].data = array
        
        if slice == SCIi:
            pixel_aligned_file[SCIi].header = copy.deepcopy(reference[SCIi].header)
        else:
            pixel_aligned_file[slice].header['NAXIS1'] = NAXIS1
            pixel_aligned_file[slice].header['NAXIS2'] = NAXIS2
            print(pixel_aligned_file[slice].header['EXTNAME'])

    pixel_aligned_file.writeto(f'{path}/{name}.fits', overwrite=True)

print(refname)
print(files[-1])
print(names[:-1])
p = Pool(7)
p.map(align, names)
#align(names[-1])