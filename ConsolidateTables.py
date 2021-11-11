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
    parser.add_argument("--list", "-l", type=str,
                        help="Path to file containing the list of files to be changed")
    parser.add_argument("--datadirectory", "-dd", type=str,
                        help="Path to folder with the data to be converted")
    parser.add_argument("--tablename", "-tn", type=str, default="LOGANTable",
                        help="Name of combined table to be generated. Without '.fits' at the end")
    parser.add_argument("--new", "-n", action='store_true',
                        help="Start new table instead of appending to table if a file with name [tablename] already exists")
    parser.add_argument("--tableindex", "-ti", type=int, default=None,
                        help="Current index of the table. If None, it will have to search the directory for the next value.")
    args = parser.parse_args()
    tablename = args.tablename
    working_directory = os.getcwd()
    if args.new:
        writingMode = "new"
    elif tablename in os.listdir(working_directory):
        writingMode = "add"
    else:
        writingMode = "other"  # TODO: Decide what is going to happen with this case

    if args.tableindex is None:
        args.tableindex = 1
        # ToDo: search directory for the correct index

    writingMode = "new"  # TODO: TESTING, REMOVE THIS
    if writingMode == "new":
        numberOfSpectraPoints = 5  # TODO: This should reflect the number of points in the spectra. For testing only 5 points will be taken
        numberOfFilters = 5  # TODO: This should reflect the number of filter. For testing only 5 points will be taken
        columnNames = ["ID", "UuidI", "UuidL", "Noise"] + ["s" + str(i) for i in range(1, numberOfSpectraPoints + 1)] + ["f" + str(i) for i in range(1, numberOfFilters + 1)]
        columnFormat = ["I", "36A", "36A", "L"] + ["E" for i in range(1, numberOfSpectraPoints + 1)] + ["E" for i in range(1, numberOfFilters + 1)]
        dummyValues = [int(1), '881ec5c0-3f57-42fa-a93b-2c70042c41k0', '1b350ffd-7df1-4d39-91ce-eqq55a48a38b', True,
                       1.368803e-17, 3.3879602e-17, 4.4652237e-17, 4.4324935e-17, 4.6498633e-17, 22.691181515985,
                       20.565539790257, 20.136398166431, 21.188410278489, 20.266081643397]
        colList = [fits.Column(name=columnNames[i], format=columnFormat[i], array=np.array([dummyValues[i]])) for i in
                   range(len(columnNames))]
        print("-------")
        coldefs = fits.ColDefs(colList)
        print(coldefs)
        dummyTable = fits.BinTableHDU.from_columns(coldefs)
    elif writingMode == "add":
        # TODO: Open file and read table
        pass

    print("before writing")
    dummyTable.writeto("tableFits2.fits", overwrite=True)
    print("after writing")

    with fits.open("tableFits2.fits") as hdul:
        data = hdul[1].data
        print(data)
        for i in data[0]:
            print(type(i))

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
    print("\n\n\n\n\n\n\n\n\n\n\n\n")
    for file in filelist:
        print(file)
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

            for k in newData:
                print(k, type(k[0]))

        file_i += 1


    print(newData)
    print("----------------------------------------------")
    columnNames = ["ID", "UuidI", "UuidL", "Noise"] + ["s" + str(i) for i in range(1, numberOfSpectraPoints + 1)] + ["f" + str(i) for i in range(1, numberOfFilters + 1)]
    columnFormat = ["I", "36A", "36A", "L"] + ["E" for i in range(1, numberOfSpectraPoints + 1)] + ["E" for i in range(1, numberOfFilters + 1)]
    colList = [fits.Column(name=columnNames[i], format=columnFormat[i], array=newData[i]) for i in range(len(columnNames))]
    coldefs = fits.ColDefs(colList)
    print(coldefs)
    dummyTable = fits.BinTableHDU.from_columns(coldefs)
    dummyTable.writeto(os.path.join(args.datadirectory, args.tablename + str(args.tableindex) + ".fits"), overwrite=True)

    # fits.append(os.path.join(args.datadirectory, args.tablename + ".fits"),)
    filename1 = os.path.join(args.datadirectory, args.tablename + ".fits")
    filename2 = os.path.join(args.datadirectory, args.tablename + str(args.tableindex) + ".fits")
    print(f"file1: {filename1}")
    print(f"file2: {filename2}")
    with fits.open(filename1) as hdul1:
        with fits.open(filename2) as hdul2:
            nrows1 = hdul1[1].data.shape[0]
            nrows2 = hdul2[1].data.shape[0]
            nrows = nrows1 + nrows2
            hdu = fits.BinTableHDU.from_columns(hdul1[1].columns, nrows=nrows)
            for colname in hdul1[1].columns.names:
                hdu.data[colname][nrows1:] = hdul2[1].data[colname]
    hdu.writeto(filename1, overwrite=True)

    # Test Reading the table appended
    print("----------------------------------------------")
    with fits.open(filename1) as hdul:
        data = hdul[1].data
        print(data.shape)
        print(data)
        print("_-_-_-_-_-_-_-_-_-_-_-_-_-_-_")
        for k in data[0]:
            print(k, type(k))



    # for file in lines:
    #     # Input
    #     print(os.path.join(args.data, "Input_" + file + ".fits"))
    #     with fits.open(os.path.join(args.data, "Input_" + file + ".fits")) as hdul:
    #         hdul.info()
    #         data = hdul[0].data.astype(np.float32)
    #         hdul[0].data = data
    #         hdul.info()
    #         # R saves scientific notation with lower case 'e', while the FITS standard indicates it should be upper case 'E'
    #         hdul.writeto(os.path.join(args.data, "Input_" + file + ".fits"), overwrite=True, output_verify="silentfix")
    #
    #     # Label
    #     print(os.path.join(args.data, "Label_" + file + ".fits"))
    #     with fits.open(os.path.join(args.data, "Label_" + file + ".fits")) as hdul:
    #         hdul.info()
    #         data = hdul[0].data.astype(np.float32)
    #         hdul[0].data = data
    #         hdul.info()
    #         # R saves scientific notation with lower case 'e', while the FITS standard indicates it should be upper case 'E'
    #         hdul.writeto(os.path.join(args.data, "Label_" + file + ".fits"), overwrite=True, output_verify="silentfix")


if __name__ == "__main__":
    consolidateFITSTables()
