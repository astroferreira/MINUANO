#!/usr/bin/env python
import numpy as np
import os
import re
import astropy.io.fits as fits
from astropy.wcs import WCS
from astropy.stats import median_absolute_deviation, biweight_location, biweight_scale
from glob import glob
from lineno import *
import inspect

#
#-------------------------------------------------------------------------------
#
def prepare_wisps(flat, median_file, output, debug):

    hdulist = fits.open(flat)
    header = hdulist[0].header
    image = hdulist['SCI'].data
    flatfield = np.fliplr(image)
    hdulist.close()
    if(debug == True):
        print("flat        ", flat)
        print("median file ", median_file)
#        
    hdulist = fits.open(median_file)
    nhdu  = len(hdulist)
    if(debug == True):
        hdulist.info()
    median_stack = hdulist['SCI'].data
    hdulist.close()

    median_ff = median_stack/flatfield

    biwt = biweight_location(median_ff,ignore_nan=True)
    print("biwt ", biwt)

    median_ff = median_ff - biwt
    fits.writeto(output,median_ff,header=header,overwrite=True)
    print("wrote file ", output)
    return
#
#===============================================================================
#

root_dir = '/home/ppxlf2/data/JWST/CEERS_REAL/'

 
rate_dir    = root_dir+'calibrated/'
skyflat_dir = root_dir+'skyflats/'
debug  = False
frame  = 1
sca = dict([
    ('481', 'NRCA1'),
    ('482', 'NRCA2'),
    ('483', 'NRCA3'),
    ('484', 'NRCA4'),
    ('485', 'NRCALONG'),
    #
    ('486', 'NRCB1'),
    ('487', 'NRCB2'),
    ('488', 'NRCB3'),
    ('489', 'NRCB4'),
    ('490', 'NRCBLONG')
    ])
flat_root = '/usr/local/nircamsuite/ncdhas_dms/cal/Flat/ISIMCV3/'
flat_suffix = 'CLEAR_2016-04-05.fits'
flat_hash = dict([
    ('NRCA1',   'NRCA1_17004_PFlat_'),
    ('NRCA2',   'NRCA2_17006_PFlat_'),
    ('NRCA3',   'NRCA3_17012_PFlat_'),
    ('NRCA4',   'NRCA4_17048_PFlat_'),
    ('NRCALONG','NRCA5_17158_PFlat_'),
    ('NRCB1',   'NRCB1_16991_PFlat_'),
    ('NRCB2',   'NRCB2_17005_PFlat_'),
    ('NRCB3',   'NRCB3_17011_PFlat_'),
    ('NRCB4',   'NRCB4_17047_PFlat_'),
    ('NRCBLONG','NRCB5_17161_PFlat_')
])
filters = []
filters.append('F150W')
filters.append('F150W2')
filters.append('F200W')
filters.append('F210M')
scas = []
scas.append(483)
scas.append(488)
scas.append(489)
wisps =[]

#for ii in range (0, len(scas)):
for ii in range (0, 1):
    detector = scas[ii]
    sca_name = sca[str(detector)]
    for jj in range(0, len(filters)):
        filter = filters[jj]
#        print(sca_name, flat_root, flat_hash[sca_name], filter, flat_suffix)
        flat   = flat_root+flat_hash[sca_name]+filter+'_'+flat_suffix
        median_file = skyflat_dir+'median_sky_'+sca_name.lower()+'_'+filter+'.fits'
        wisp_file = re.sub('median_sky','wisps',median_file)
        print("flat       is ", flat)
        print("median sky is ", median_file)
        print("wisp       is ", wisp_file)
        prepare_wisps(flat, median_file, wisp_file, debug)

files = sorted(glob(rate_dir+'*nrcb4_rate.fits'))
#print(files)
for index in range(len(files)):
    file = files[index]
    hdulist = fits.open(file)
    nhdu  = len(hdulist)
    if(debug == 1):
        hdulist.info()

    for ii in range(0, nhdu):
        header = hdulist[ii].header
        if('DETECTOR' in header):
            detector = header['detector']
        if('FILTER' in header):
            filter = header['FILTER']

    print("filter ", filter, 'detector ', detector)
    org =hdulist['SCI'].data
    hdulist.close()
    #
# The naming is inconsistent between rate files and flatfields
#
    sca_name = detector
    wisp_file = skyflat_dir+'wisps_'+sca_name.lower()+'_'+filter+'.fits'
#
    wisplist   = fits.open(wisp_file)
    wisplist.info()
    median_ff =  wisplist[0].data
    wisplist.close()

    corrected = org - median_ff
    new_file = re.sub('rate.fits','corr.fits',file)
    print("org file", file)
    print("new file", new_file)
    fits.writeto(new_file,corrected, header, overwrite=True)
#    exit(0)
    
