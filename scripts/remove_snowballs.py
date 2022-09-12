from astropy.io import fits
import sys


def mask_snowballs(file, snowball_erode=3, snowball_dilate=18, mask_bit=1024, instruments=['NIRCAM'], unset4=False, **kwargs):

    import scipy.ndimage as nd
    
    
    with fits.open(file, mode='update') as _im:
    
        crs = (_im['DQ'].data & 4)
        
        _erode = nd.binary_erosion(crs > 0, iterations=snowball_erode)
        snowball_mask = nd.binary_dilation(_erode, 
                                            iterations=snowball_dilate)
        
        label, num_labels = nd.label(snowball_mask)
        
        _im['DQ'].data |= (snowball_mask*mask_bit).astype(_im['DQ'].data.dtype)
        
        if unset4:
            _im['DQ'].data -= (_im['DQ'].data & 4)
            
        _im['SCI'].header['SNOWMASK'] = (True, 'Snowball mask applied')
        _im['SCI'].header['SNOWEROD'] = (snowball_erode,
                                            'CR bit=4 erosion for snowballs')
        _im['SCI'].header['SNOWDILA'] = (snowball_dilate,
                                        'CR bit=4 dilation for snowballs')
        _im['SCI'].header['SNOWMBIT'] = (mask_bit, 'Snowball bit ')
        _im['SCI'].header['SNOWBALN'] = (num_labels,
                            'Number of labeled features in snowball mask')
        
        _maskfrac = snowball_mask.sum() / snowball_mask.size
        _im['SCI'].header['SNOWBALF'] = (_maskfrac,
                                'Fraction of masked pixels in snowball mask')
        
        msg = f"Snowball mask: {file} "
        msg += f" N={_im['SCI'].header['SNOWBALN']:>3}"
        msg += f"(f={_im['SCI'].header['SNOWBALF']*100:.2f}%)"
        


        _im.flush()

    return True

file = sys.argv[1]
mask_snowballs(file)