import argparse
import os
from astropy.io import fits
import numpy as np
import pandas as pd


def consolidate_files_into_tables():
    """
    Generates a single FITS table containing the information of the whole set of instructions
    """
    # Initialize parameters, read argparse and reads files

    parser = argparse.ArgumentParser()
    parser.add_argument("--list", "-l", type=str, default="ListToConvertoTo32.txt",
                        help="Path to file containing the list of files to be changed")
    parser.add_argument("--datadirectory", "-dd", type=str,
                        help="Path to folder with the data to be converted")
    parser.add_argument("--nametable", "-nt", type=str, default="LOGANTable.fits",
                        help="Name of combined table to be generated. If absolute path is NOT given,"
                             "will search for it in the data directory>current working directory."
                             "If the file is not found, it will be generated in the Data directory")
    parser.add_argument("--namelabel", "-nl", type=str, default="LOGANLables.fits",
                        help="Name of combined lables to be generated. If absolute path is NOT given,"
                             "will search for it in the data directory>current working directory."
                             "If the file is not found, it will be generated in the Data directory")
    parser.add_argument("--namemetadata", "-nm", type=str, default="LOGANMetaData.fits",
                        help="Name of combined metadata to be generated. If absolute path is NOT given,"
                             "will search for it in the data directory>current working directory."
                             "If the file is not found, it will be generated in the Data directory")
    parser.add_argument("--new", "-n", action='store_true',
                        help="Start new table instead of appending to table if a file with name"
                             "[nametable] already exists")
    parser.add_argument("-ns", type=int, help="Number of points in the spectra.")
    parser.add_argument("-nf", type=int, help="Number of filters/photometry measured.")
    parser.add_argument("-nid", type=int, help="ID of last element.")
    args = parser.parse_args()
    # TODO: Do the equivalent for Labels and metadata!

    # Locate the Consolidated table
    # TODO: Optional: do the same for list
    list_filename, _ = search_file(args, args.list, "File list")
    name_table, writing_table = search_file(args, args.nametable, "Table")
    # name_label, writing_lable = search_file(args, args.namelabel, "Labels")
    # name_metadata, writing_metadata = search_file(args, args.namemetadata, "Metadata")

    if args.new:
        print(f"Writing mode forced to generate a new file.")
        writing_table = "new"

    column_names = ["ID", "UuidI", "UuidL", "Noise"] +\
                   ["s" + str(i) for i in range(1, args.ns + 1)] +\
                   ["f" + str(i) for i in range(1, args.nf + 1)]
    column_format = ["I", "36A", "36A", "L"] +\
                    ["E" for _ in range(1, args.ns + 1)] +\
                    ["E" for _ in range(1, args.nf + 1)]
    if writing_table == "new":
        # If new, generate a dummy (empty) table.
        col_list = [fits.Column(name=column_names[i], format=column_format[i]) for i in range(len(column_names))]
        col_defs = fits.ColDefs(col_list)
        tmp_table = fits.BinTableHDU.from_columns(col_defs)
        tmp_table.writeto(name_table, overwrite=True)

    # Read the file with the names of files to be read and generate filelist
    with open(list_filename) as file:
        filelist = file.readlines()
        filelist = [line.rstrip().rstrip("\x00") for line in filelist]

    number_of_rows = len(filelist)
    new_data = [np.zeros(number_of_rows, dtype=np.int16),
                np.zeros(number_of_rows, dtype=np.dtype('U36')),
                np.zeros(number_of_rows, dtype=np.dtype('U36')),
                np.zeros(number_of_rows, dtype=np.bool_)] + \
               [np.zeros(number_of_rows, dtype=np.float32)] * args.ns + \
               [np.zeros(number_of_rows, dtype=np.float32)] * args.nf

    file_i = 0
    last_id = args.nid
    for file in filelist:
        with fits.open(os.path.join(args.datadirectory, "Input_" + file + ".fits")) as hdul:
            data_file = hdul[0].data
            hdr = hdul[0].header

            # ===========================================================
            # Insert new information into its place
            # ID
            new_data[0][file_i] = last_id + 1
            last_id += 1
            # UUIDs
            new_data[1][file_i] = hdr["UUIDINP"]
            new_data[2][file_i] = hdr["UUIDLAB"]
            # Noise
            new_data[3][file_i] = hdr["Noise"]
            # Spectra
            for s in range(4, 4 + args.ns):
                new_data[s][file_i] = data_file[0][s - 4]
            # Filters
            for f in range(4 + args.ns, 4 + args.ns + args.nf):
                new_data[f][file_i] = hdr["FILTERV" + str(f - 3 - args.ns)]
        file_i += 1

    # Transform the matrix into a temporary file
    col_list = [fits.Column(name=column_names[i], format=column_format[i], array=new_data[i])
                for i in range(len(column_names))]
    col_defs = fits.ColDefs(col_list)
    tmp_table = fits.BinTableHDU.from_columns(col_defs)
    tmpname = os.path.join(args.datadirectory, "tmp.fits")
    tmp_table.writeto(tmpname, overwrite=True)

    with fits.open(name_table) as hdul1:
        with fits.open(tmpname) as hdul2:
            nrows1 = hdul1[1].data.shape[0]
            nrows2 = hdul2[1].data.shape[0]
            nrows = nrows1 + nrows2
            print(nrows1, nrows2, nrows)
            hdu = fits.BinTableHDU.from_columns(hdul1[1].columns, nrows=nrows)
            for colname in hdul1[1].columns.names:
                hdu.data[colname][nrows1:] = hdul2[1].data[colname]
    hdu.writeto(name_table, overwrite=True)

    # Remove files that have been read and stored (tmp, list-file, and read files)
    if False:       # ToDo: Made inaccesible during testing
        os.remove(os.path.join(args.datadirectory, "tmp.fits"))
        for file in filelist:
            os.remove(os.path.join(args.datadirectory, "Input_" + file + ".fits"))
            os.remove(os.path.join(args.datadirectory, "Label_" + file + ".fits"))
        os.remove(list_filename)

    # TODO: Remove tmp file and fits files

    # Test Reading the table appended
    # TODO: This is temporarily stored here for future use when reading the table
    # TODO: For future usage, store in NPY File
    # TODO: https://machinelearningmastery.com/how-to-save-a-numpy-array-to-file-for-machine-learning/
    # print("----------------------------------------------")
    # with fits.open(name_table) as hdul:
    #     data = hdul[1].data
    #     print(data.shape)
    #     print("_-_-_-_-_-_-_-_-_-_-_-_-_-_-_")
    #     for k in data[0]:
    #         print(k, type(k))
    #
    #     print("@@@@@@@@@@@")
    #     print(data.shape)
    #     names = hdul[1].columns.names  # we need the column names
    #     cols = [hdul[1].data.field(col) for col in names]  # and their content
    #     cat = np.rec.fromarrays(cols, names=names)
    #     print(cat)
    #     print(cat.dtype.names)  # similar to hdus[1].columns.names
    #     print(cat.shape)  # and we have access to numpy commands
    #
    #     selection_names = ["s" + str(i) for i in range(1, args.ns+1)] + ["f" + str(i) for i in range(1, args.nf + 1)]
    #     selection = cat[selection_names]
    #     print(selection)
    #     print(selection.shape)
    #     print(type(selection))
    #     arr = pd.DataFrame(selection).to_numpy()
    #     print(arr)
    #     print(arr.shape)
    #     print(type(arr))
    #     print(type(arr[0][0]))


def search_file(args, tablename, desc):
    if not os.path.isabs(tablename):
        # If the name is NOT an absolute path, search in the data directory and then in the cwd
        writing_mode = "add"
        if tablename in os.listdir(args.datadirectory):
            tablename = os.path.join(args.datadirectory, tablename)
            print(f"{desc} was found in {tablename}")
        elif tablename in os.listdir(os.getcwd()):
            tablename = os.path.join(os.getcwd(), tablename)
            print(f"{desc} was found in {tablename}")
        else:
            tablename = os.path.join(args.datadirectory, tablename)
            writing_mode = "new"
            print(f"{desc} was not found. It will be generated in {tablename}")
    else:
        if os.path.isfile(tablename):
            print(f"{desc} was found in {tablename}")
            writing_mode = "add"
        else:
            print(f"{desc} was not found. It will be generated in {tablename}")
            writing_mode = "new"
    return tablename, writing_mode


if __name__ == "__main__":
    print("\n\n\n")
    consolidate_files_into_tables()
