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
    filename = ls[0]

    with fits.open(os.path.join(filedir, filename)) as hdul:
        hdul.info()
        data = hdul[0].data
        print(type(data))
        data = np.float32(data)
        print(type(data))
        hdul[0].data = data
        hdul.info()
        hdul.writeto(os.path.join(filedir, "aaaaa.fits"), overwrite=True)





if __name__ == "__main__":
    main()
