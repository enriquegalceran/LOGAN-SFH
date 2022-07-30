# -*- coding: utf-8 -*-
# SpeCtra Observer and TesTer [Cyclops]

# import json
import os
# import typing
# import subprocess
# import uuid
import argparse
# from shutil import copyfile
# from datetime import datetime
# import numpy as np
# import keras
import matplotlib.pyplot as plt
from astropy.io import fits
import time
import matplotlib
print(matplotlib.matplotlib_fname())

def main(update_values=None, **mainkwargs):
    # Argument Parser (argparse)
    parser = argparse.ArgumentParser(description="Spectra and Data Observer")

    parser.add_argument("-ln", "--logname", default=None, type=str,
                        help="If a filename is given, it will clean this log (i.e. remove extra characters")
    parser.add_argument("-clo", "--cleanlogoutput", default=None, type=str,
                        help="If a filename is given, will write the cleaned log in this location. "
                             "Requires a value in --logname to work.")
    parser.add_argument("-vm", "--verify_model", default=None, type=list, nargs="*",
                        help="Test the prediction of the Neural Network. the absolute/relative path to the model is"
                             " required. Requieres a dataset to verify (--data_path).")
    parser.add_argument("-dp", "--data_path", default=None, type=str,
                        help="Indicates a dataset to verify the model. If the data is not a path (absolute/relative) "
                             "(if it does not end in '.fits'), it will consider the data to be located in the current "
                             "folder and that the given value is the sufix. REQUIRES THE INPUT DATA. The Labels are "
                             "expected to be located in the same folder.")
    parser.add_argument("-v", "--verbose", default=1, type=int,
                        help="Indicates the value of verbosity. Minimum=0, Default=1.")
    parser.add_argument("-jc", "--json_clean", action="append", nargs="+", default=None, type=str,
                        help="Path to json file to be cleaned. Multiple files can be given (either using -jc again, or "
                             "be simply giving multiple paths).If verbose > 0, json file will be printed.")
    args = parser.parse_args()

    with fits.open("/Users/enrique/Documents/GitHub/LOGAN-SFH/Testing_spectra1.fits") as hdul1:
        data1 = hdul1[0].data
    with fits.open("/Users/enrique/Documents/GitHub/LOGAN-SFH/Testing_spectra2.fits") as hdul1:
        data2 = hdul1[0].data


    color = "orange"
    plt.ion()
    fig, ax = plt.subplots()
    ln, = ax.plot(data1)
    while True:
        ln.set_color(color)
        # plt.show(block=False)
        print("asdfg")
        time.sleep(1)
        color = input("Color a cambiar")
        if color == "end":
            break


if __name__ == "__main__":
    # with fits.open("/Users/enrique/Documents/GitHub/LOGAN-SFH/TrainingData/Input_combined_wave_steps.fits") as hdul:
    #     spectra = hdul[0].data
    #     output1 = spectra[1:-4, 1]
    #     output2 = spectra[1:-4, 2]
    #
    # hdu = fits.PrimaryHDU(output1)
    # hdu.writeto("/Users/enrique/Documents/GitHub/LOGAN-SFH/Testing_spectra1.fits")
    # hdu = fits.PrimaryHDU(output2)
    # hdu.writeto("/Users/enrique/Documents/GitHub/LOGAN-SFH/Testing_spectra2.fits")
    main()
