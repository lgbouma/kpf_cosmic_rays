from astropy.io import fits

def read_fits(fits_file,ext=0):
    '''
    Shortcut function to get the header and data from a fits file and a given
    extension.
    '''

    hdulist = fits.open(fits_file)
    img_header = hdulist[ext].header
    img_data = hdulist[ext].data
    hdulist.close()

    return img_data, img_header


