import numpy as np
import os
import re
import sys
import astropy.io.fits as fits
from os.path import exists
from glob import glob
import inspect
#
#==============================================================================
# Wisps are seen in F150W, F150W2, F200W, F210M
# Status unknown in F162M
# F182M needs a scaling factor
# SCAs presenting significant wisps : A3, A4, B3, B4
template_dir  = '/home/ppxlf2/data/JWST/MINUANO/wisp_templates/'
template = dict([
    ('NRCA3_F150W', 'wisps_nrca3_F150W.fits'),
    ('NRCA3_F150W2','wisps_nrca3_F150W2.fits'),
    ('NRCA3_F182M', 'wisps_nrca3_F200W.fits'),
    ('NRCA3_F200W', 'wisps_nrca3_F200W.fits'),
    ('NRCA3_F210M', 'wisps_nrca3_F210M.fits'),
    ('NRCA4_F150W', 'wisps_nrca4_F150W.fits'),
    ('NRCA4_F150W2','wisps_nrca4_F150W2.fits'),
    ('NRCA4_F182M', 'wisps_nrca4_F200W.fits'),
    ('NRCA4_F200W', 'wisps_nrca4_F200W.fits'),
    ('NRCA4_F210M', 'wisps_nrca4_F210M.fits'),
    ('NRCB3_F150W', 'wisps_nrcb3_F150W.fits'),
    ('NRCB3_F150W2','wisps_nrcb3_F150W2.fits'),
    ('NRCB3_F182M', 'wisps_nrcb3_F200W.fits'),
    ('NRCB3_F200W', 'wisps_nrcb3_F200W.fits'),
    ('NRCB3_F210M', 'wisps_nrcb3_F210M.fits'),
    ('NRCB4_F150W', 'wisps_nrcb4_F150W.fits'),
    ('NRCB4_F150W2','wisps_nrcb4_F150W2.fits'),
    ('NRCB4_F182M', 'wisps_nrcb4_F210M.fits'),
    ('NRCB4_F200W', 'wisps_nrcb4_F200W.fits'),
    ('NRCB4_F210M', 'wisps_nrcb4_F210M.fits')
    ])
    
#
debug  = 0
overwrite = False
overwrite = True
rate_dir = sys.argv[1]
#print("rate_dir is ", rate_dir)
#files = sorted(glob(rate_dir+'*rate.fits'))
#print(files)
#for index in range(len(files)):
file = rate_dir#files[index]
hdulist = fits.open(file)
nhdu  = len(hdulist)
wispfile = None
if(debug == 1):
    hdulist.info()

for ii in range(0, nhdu):
    header = hdulist[ii].header
    if('DETECTOR' in header):
        detector = header['detector']
    if('FILTER' in header):
        filter = header['FILTER']
    if('WISPSUB' in header):
        wispfile = header['WISPSUB']
#    hdulist.close()
        

#
# The very few tests carried out in commissioning showed that a
# plain subtraction was effective. If this is not the case, it may
# be necessary toscale the template intensity to determine 
# an optimal subtractionlevel
#os.path.exists(path_to_file)
print("filter ", filter, 'detector ', detector)
if(filter == 'F150W' or filter == 'F150W2' or \
    filter == 'F182M' or \
    filter == 'F200W' or filter == 'F210M'):
    if(wispfile != None):
        exit()
        hdulist.close()
        corr_file = re.sub('rate.fits','rate_corr.fits',file)
        command = 'ln -s '+file +' '+corr_file
        print(command)
        os.system(command)
        #continue
    if(filter == 'F182M'):
        scale = 1.50
    else:
        scale = 1.00
    data = np.array(hdulist[1].data)
    print("file ", file)
    corr_file = re.sub('rate.fits','rate_corr.fits',file)
    if (os.path.exists(corr_file) and overwrite == False):
        print('corrected file exists and being skipped ', file)
        exit()#continue
    key = detector+'_'+filter
    try:
        wisp = template_dir+template[key]
    except:
        exit()#continue
    print("wisp ", wisp)
    median = fits.getdata(wisp,ext=0)
    #median = median[::-1, ::-1].copy()
    # Should one mask or massage the template NaNs ?
    yy, xx = np.where(np.isnan(median))
    # 
    # mask   = np.zeros(median.shape)
    # mask[yy,xx] = 1
    # new_mask = np.where(mask > 0, True, False)
    # mask = np.array(new_mask)
    # corrected = np.ma.array(data) - np.ma.array(median,mask=mask)*scale
    # corrected = np.array(corrected)
    # we went for massaging...
    massaged = np.array(median)
    massaged[yy,xx] = 0
    corrected = data - massaged*scale
    path = wisp.split('/')
    wispfile = path[len(path)-1]
    hdulist[0].header['WISPSUB'] = wispfile
    hdulist[1].data = corrected
    print("corrected file is ", corr_file)
    hdulist.writeto(corr_file,overwrite=overwrite)
    hdulist.close()
#        exit(0)