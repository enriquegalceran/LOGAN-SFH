import sys
import argparse
import os
from astropy.io import fits
import numpy as np


def convert_bytes(num):
    """
    this function will convert bytes to MB.... GB... etc
    """
    for x in ['bytes', 'KB', 'MB', 'GB', 'TB', 'PB']:
        if num < 1024.0:
            return "%3.1f %s" % (num, x)
        num /= 1024.0

def file_size(file_path, convert=True):
    """
    this function will return the file size
    """
    if os.path.isfile(file_path):
        file_info = os.stat(file_path)
        if convert:
            return convert_bytes(file_info.st_size)
        else:
            return file_info.st_size


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-list", type=str, help="Path to file containing the list of files to be changed")
    parser.add_argument("-data", type=str, help="Path to folder with the data to be converted")
    args = parser.parse_args()
    print(args)
    filename = args.list
    print(filename)
    with open(filename) as file:
        lines = file.readlines()
        lines = [line.rstrip().rstrip("\x00") for line in lines]

    print(f"Total amount of files to be converted: {len(lines)}")

    for file in lines:
        # Input
        print(os.path.join(args.data, "Input_" + file + ".fits"))
        print(f"Initial filesize: {file_size(os.path.join(args.data, 'Input_' + file + '.fits'), False)}")
        with fits.open(os.path.join(args.data, "Input_" + file + ".fits")) as hdul:
            hdul.info()
            data = hdul[0].data.astype(np.float32)
            hdul[0].data = data
            hdul.info()
            # R saves scientific notation with lower case 'e', while the FITS standard indicates it should be upper case 'E'
            hdul.writeto(os.path.join(args.data, "Input_" + file + ".fits"), overwrite=True, output_verify="silentfix")
        print(f"Final filesize: {file_size(os.path.join(args.data, 'Input_' + file + '.fits'), False)}")

        # Label
        print(os.path.join(args.data, "Label_" + file + ".fits"))
        print(f"Initial filesize: {file_size(os.path.join(args.data, 'Label_' + file + '.fits'), False)}")
        with fits.open(os.path.join(args.data, "Label_" + file + ".fits")) as hdul:
            hdul.info()
            data = hdul[0].data.astype(np.float32)
            hdul[0].data = data
            hdul.info()
            # R saves scientific notation with lower case 'e', while the FITS standard indicates it should be upper case 'E'
            hdul.writeto(os.path.join(args.data, "Label_" + file + ".fits"), overwrite=True, output_verify="silentfix")
        print(f"Final filesize: {file_size(os.path.join(args.data, 'Label_' + file + '.fits'), False)}")

    # Remove file
    os.remove(filename)





if __name__ == "__main__":
    main()
