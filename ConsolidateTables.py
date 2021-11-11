import sys
import argparse
import os
from astropy.io import fits
from astropy.table import Table, vstack
import numpy as np


def consolidateFITSTables():
    """
    Generates a single FITS table containing the information of the whole set of instructions
    """
    # Initialize parameters, read argparse and reads files

    parser = argparse.ArgumentParser()
    parser.add_argument("--list", "-l", type=str, default="ListToConvertoTo32.txt",
                        help="Path to file containing the list of files to be changed")
    parser.add_argument("--datadirectory", "-dd", type=str,
                        help="Path to folder with the data to be converted")
    parser.add_argument("--tablename", "-tn", type=str, default="LOGANTable.fits",
                        help="Name of combined table to be generated. If absolute path is NOT given,"
                             "will search for it in the data directory>current working directory."
                             "If the file is not found, it will be generated in the Data directory")
    parser.add_argument("--new", "-n", action='store_true',
                        help="Start new table instead of appending to table if a file with name [tablename] already exists")
    args = parser.parse_args()

    # Locate the Consolidated table
    tablename = args.tablename
    if not os.path.isabs(tablename):
        # If the name is NOT an absolute path, search in the data directory and then in the cwd
        writingMode = "add"
        if tablename in os.listdir(args.datadirectory):
            tablename = os.path.join(args.datadirectory, tablename)
            print(f"Table was found in {tablename}")
        elif tablename in os.listdir(os.getcwd()):
            tablename = os.path.join(os.getcwd(), tablename)
            print(f"Table was found in {tablename}")
        else:
            tablename = os.path.join(args.datadirectory, tablename)
            writingMode = "new"
            print(f"Table was not found. It will be generated in {tablename}")
    else:
        if os.path.isfile(tablename):
            print(f"Table was found in {tablename}")
            writingMode = "add"
        else:
            print(f"Table was not found. It will be generated in {tablename}")
            writingMode = "new"

    if args.new:
        print(f"Writing mode forced to generate a new file.")
        writingMode = "new"

    numberOfSpectraPoints = 5  # TODO: This should reflect the number of points in the spectra. For testing only 5 points will be taken
    numberOfFilters = 5  # TODO: This should reflect the number of filter. For testing only 5 points will be taken
    columnNames = ["ID", "UuidI", "UuidL", "Noise"] + ["s" + str(i) for i in range(1, numberOfSpectraPoints + 1)] + ["f" + str(i) for i in range(1, numberOfFilters + 1)]
    columnFormat = ["I", "36A", "36A", "L"] + ["E" for i in range(1, numberOfSpectraPoints + 1)] + ["E" for i in range(1, numberOfFilters + 1)]
    if writingMode == "new":
        # If new, generate a dummy (empty) table.
        colList = [fits.Column(name=columnNames[i], format=columnFormat[i]) for i in range(len(columnNames))]
        coldefs = fits.ColDefs(colList)
        dummyTable = fits.BinTableHDU.from_columns(coldefs)
        dummyTable.writeto(tablename, overwrite=True)
    elif writingMode == "add":
        # TODO: Open file and read table
        pass

    # Read the file with the names of files to be read and generate filelist
    listfilename = args.list
    with open(listfilename) as file:
        filelist = file.readlines()
        filelist = [line.rstrip().rstrip("\x00") for line in filelist]

    filelist = ["testFits1.fits", "testFits2.fits", "testFits3.fits", "testFits4.fits", "testFits5.fits", "testFits6.fits"]

    # numberOfRows = len(filelist)  # ToDo: This will be the way to calculate numberOfRows
    numberOfRows = 5
    newData = [np.zeros(numberOfRows, dtype=np.int16),
               np.zeros(numberOfRows, dtype=np.dtype('U36')),
               np.zeros(numberOfRows, dtype=np.dtype('U36')),
               np.zeros(numberOfRows, dtype=np.bool_)] + \
              [np.zeros(numberOfRows, dtype=np.float32)] * numberOfSpectraPoints + \
              [np.zeros(numberOfRows, dtype=np.float32)] * numberOfFilters

    file_i = 0
    lastID = 10  # TODO: This needs to be read from the previous existing table or start at 1 (0?)
    # print("\n\n\n\n\n\n\n\n\n\n\n\n")
    for file in filelist:
        # print(file)
        if file == filelist[numberOfRows]:  # TODO: Just for testing, this will be removed
            break

        # with fits.open(os.path.join(args.datadirectory, "Input_" + file + ".fits")) as hdul:
        with fits.open(os.path.join(args.datadirectory, file)) as hdul:
            dataFile = hdul[0].data
            hdr = hdul[0].header

            # ===========================================================
            # Insert new information into its place
            # ID
            newData[0][file_i] = lastID + 1
            lastID += 1
            # UUIDs
            newData[1][file_i] = hdr["UUIDINP"]
            newData[2][file_i] = hdr["UUIDLAB"]
            # Noise
            newData[3][file_i] = hdr["Noise"]
            # Spectra
            # TODO: This is truncated for testing. needs to be increased accordingly
            for s in range(4, 4 + numberOfSpectraPoints):
                newData[s][file_i] = dataFile[0][s - 4]
            # Filters
            for f in range(4 + numberOfSpectraPoints, 4 + numberOfSpectraPoints + numberOfFilters):
                newData[f][file_i] = hdr["FILTERV" + str(f - 3 - numberOfSpectraPoints)]

            # for k in newData:
            #     print(k, type(k[0]))

        file_i += 1


    colList = [fits.Column(name=columnNames[i], format=columnFormat[i], array=newData[i]) for i in range(len(columnNames))]
    coldefs = fits.ColDefs(colList)
    dummyTable = fits.BinTableHDU.from_columns(coldefs)
    dummyTable.writeto(os.path.join(args.datadirectory, "tmp.fits"), overwrite=True)

    # fits.append(os.path.join(args.datadirectory, args.tablename),)
    tmpname = os.path.join(args.datadirectory, "tmp.fits")
    print(f"file1: {tablename}")
    print(f"file2: {tmpname}")
    with fits.open(tablename) as hdul1:
        with fits.open(tmpname) as hdul2:
            nrows1 = hdul1[1].data.shape[0]
            nrows2 = hdul2[1].data.shape[0]
            nrows = nrows1 + nrows2
            print(nrows1, nrows2, nrows)
            hdu = fits.BinTableHDU.from_columns(hdul1[1].columns, nrows=nrows)
            for colname in hdul1[1].columns.names:
                hdu.data[colname][nrows1:] = hdul2[1].data[colname]
    hdu.writeto(tablename, overwrite=True)

    # TODO: Remove tmp file and fits files

    # Test Reading the table appended
    print("----------------------------------------------")
    with fits.open(tablename) as hdul:
        data = hdul[1].data
        print(data.shape)
        print("_-_-_-_-_-_-_-_-_-_-_-_-_-_-_")
        for k in data[0]:
            print(k, type(k))


if __name__ == "__main__":
    print("\n\n\n")
    consolidateFITSTables()
