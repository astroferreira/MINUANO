from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from matplotlib import pyplot as plt

from astropy.visualization.mpl_normalize import simple_norm,ImageNormalize
from astropy.visualization import SqrtStretch
from astropy.visualization import simple_norm

from scipy.ndimage import zoom
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
import logging 
import argparse
import glob
import os
import numpy as np


FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT, level=logging.DEBUG)
#logging.basicConfig(filename='example.log', encoding='utf-8', )

parser = argparse.ArgumentParser(description='Create individual exposures stamps from RA/DEC')
parser.add_argument('FIELD', type=str)
parser.add_argument('RA', type=float)
parser.add_argument('DEC', type=float)
parser.add_argument('ID', type=str)

norm = ImageNormalize(stretch=SqrtStretch(), vmin=-0.00, vmax=0.5)


if __name__ == '__main__':

    args = parser.parse_args()
    RA = args.RA
    DEC = args.DEC
    FIELD = args.FIELD
    ID = args.ID

    base_path = f'/home/ppxlf2/data/JWST/{FIELD}/calibrated/'
    coord = SkyCoord(ra=RA, dec=DEC, unit=u.deg)
    
    logging.info(f'Looking source(RA={RA}, DEC={DEC}) on the {FIELD} field')

    cal_files = sorted(glob.glob(os.path.join(base_path, '*skymatchstep.fits')))
    logging.info(f'{len(cal_files)} cal files found')

    fig = plt.figure( facecolor='white', dpi=200)
    ncols = 5
    gs = fig.add_gridspec(ncols, ncols, wspace=0, hspace=0.5)
    #gs_sub = gs[0].subgridspec(10, 7, hspace=0, wspace=-0.2)

    fig.suptitle(f'ID:{ID} - {FIELD} - RA={RA}, DEC={DEC}')

    present = []
    index = 0
    cutouts = []
    stacked = np.zeros((64, 64))
    for cal in cal_files:
        data = fits.open(cal)
        wcs = WCS(data['SCI'].header)

        
        
        if coord.contained_by(wcs):
            xdx = int(index / ncols)
            ydx = int(index % ncols)

            print(index, xdx, ydx)
            
            ax = fig.add_subplot(gs[xdx, ydx])

            filter = data[0].header['FILTER']
            #module = data[0].header['MODULE']
            detector = data[0].header['DETECTOR']

            size = 32 if 'LONG' in detector else 64
            
            logging.info(f'Extracting {size}x{size} {filter}-{detector}')

            cutout = Cutout2D(data[1].data, coord, size=size, wcs=wcs)
            present.append(cal)
            ax.set_title(f'{filter}-{detector}', fontsize=6)
            ax.imshow(cutout.data, norm=norm, cmap='gist_gray_r')
            ax.set(xticks=[], yticks=[])
            ax.grid('off')
            
            if size==64:
                stacked += cutout.data
            else:
                print(size)
                flux = cutout.data.sum()
                temp = zoom(cutout.data, 2)
                temp /= stacked.sum()
                temp *= flux
                stacked += temp

            index += 1
        
    xdx = int(index / ncols)
    ydx = int(index % ncols)

            
    ax = fig.add_subplot(gs[xdx, ydx])
    ax.set_title(f'STACKED', fontsize=6)
    ax.imshow(stacked/index, cmap='gist_gray_r')
    ax.set(xticks=[], yticks=[])
    ax.grid('off')



    plt.savefig(f'{ID}.png')

           
                        
    logging.info(f'Source found in {len(present)} files')

    






    #print(cal_files)




