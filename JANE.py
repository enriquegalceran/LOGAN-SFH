# Just A Network Executor

import numpy as np
import pandas as pd
import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
import os


def loadfiles():
    # ToDo: argparse these variables
    input_path = "/Volumes/Elements/Outputs/Input_snorm_20211203T112003_mwIBuq.fits"
    labels_path = "/Volumes/Elements/Outputs/Label_snorm_20211203T112003_mwIBuq.fits"

    # Load files
    print("Loading input...")
    with fits.open(input_path) as hdul:
        input_data = hdul[0].data
        input_header = hdul[0].header

    print("Loading labels...")
    with fits.open(labels_path) as hdul:
        labels_data = hdul[0].data
        labels_header = hdul[0].header

    print(input_data.shape)


def main():
    loadfiles()


if __name__ == "__main__":
    main()
