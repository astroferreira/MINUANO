from image1overf import sub1fimaging
from astropy.io import fits


import sys

cal2file = sys.argv[1]

cal21overffile = cal2file.replace('_cal.fits','_cal_1overf.fits')
print ('Running 1/f correction on {} to produce {}'.format(cal2file,cal21overffile))
with fits.open(cal2file) as cal2hdulist:
    if cal2hdulist['PRIMARY'].header['SUBARRAY']=='FULL' or cal2hdulist['PRIMARY'].header['SUBARRAY']=='SUB256':
        sigma_bgmask=3.0
        sigma_1fmask=2.0
        splitamps=True   #Set to True only in a sparse field so each amplifier will be fit separately. 
        slowaxis = abs(cal2hdulist['PRIMARY'].header['SLOWAXIS'])
        correcteddata = sub1fimaging(cal2hdulist,sigma_bgmask,sigma_1fmask,splitamps, slowaxis=slowaxis)
        if cal2hdulist['PRIMARY'].header['SUBARRAY']=='FULL':
            cal2hdulist['SCI'].data[4:2044,4:2044] = correcteddata  
        elif cal2hdulist['PRIMARY'].header['SUBARRAY']=='SUB256':
            cal2hdulist['SCI'].data[:252,:252] = correcteddata
        """
        print(slowaxis)
        if slowaxis==1:
            slowaxis=2
        elif slowaxis==2:
            slowaxis=1

        print(slowaxis)
        correcteddata = sub1fimaging(cal2hdulist,sigma_bgmask,sigma_1fmask,splitamps, slowaxis=slowaxis)
        if cal2hdulist['PRIMARY'].header['SUBARRAY']=='FULL':
            cal2hdulist['SCI'].data[4:2044,4:2044] = correcteddata  
        elif cal2hdulist['PRIMARY'].header['SUBARRAY']=='SUB256':
            cal2hdulist['SCI'].data[:252,:252] = correcteddata
        """
        cal2hdulist.writeto(cal21overffile, overwrite=True)