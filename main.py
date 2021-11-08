import astropy
from astropy.io import fits
import numpy as np
import pandas as pd
import rpy2
import sys
import os


def main():
    print("Hello World!")

    # Open file
    filedir = "/Volumes/Elements/Outputs"
    ls = os.listdir(filedir)
    print(ls)
    filename = ls[2]

    with fits.open(os.path.join(filedir, filename)) as hdul:
        hdul.info()
        data = hdul[0].data.astype(np.float32)
        hdul[0].data = data
        # R saves scientific notation with lower case 'e', while the FITS standard indicates it should be upper case 'E'
        hdul.writeto(os.path.join(filedir, "aaaaa.fits"), overwrite=True, output_verify="silentfix")





if __name__ == "__main__":
    main()
