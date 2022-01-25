# -*- coding: utf-8 -*-
# TESts, Scripts and Auxilliary functions [Sage]

import math


def convert_bytes(size_bytes, decimals: int = 2, base_mult: int = 1024):
    """
    Converts bytes (B) to KB, MB, GB, ... up to YB
    :param size_bytes:
    :param decimals:
    :param base_mult:
    :return:
    """
    if size_bytes == 0:
        return "0B"
    size_name = ['B', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB', 'ZB', 'YB']
    i = int(math.floor(math.log(size_bytes, base_mult)))
    p = math.pow(base_mult, i)
    s = round(size_bytes / p, decimals)
    return "%s %s" % (s, size_name[i])


def print_train_test_val_sizes(split1, split2, main_title=""):
    """
    Outputs a table with the sizes of the different groups. Useful for debug.
    Inputs are first split (train+val, test) and second split (train, val).
    :param split1:
    :param split2:
    :param main_title:
    :return:
    """
    data = [
        [' ', 'train', 'val', 'test'],
        ['Spect', str(split2[0].shape), str(split2[1].shape), str(split1[1].shape)],
        ['Magnitudes', str(split2[2].shape), str(split2[3].shape), str(split1[3].shape)],
        ['SFH', str(split2[4].shape), str(split2[5].shape), str(split1[5].shape)],
        ['Metallicity', str(split2[6].shape), str(split2[7].shape), str(split1[7].shape)]
    ]
    print_table(data, main_title)


def print_table(tabl_data, main_title='', col_separation='| ', min_length=0):
    """
    Given a table (list of lists), prints a table
    :param tabl_data:
    :param main_title:
    :param col_separation:
    :param min_length:
    :return:
    """
    class Format:
        end = '\033[0m'
        underline = '\033[4m'
    ncolumns = len(tabl_data[0])
    for row in tabl_data:
        assert len(row) == ncolumns, "All rows need to have the same length."
    char_in_column = max(max([len(element) for row in tabl_data for element in row]), min_length)

    string_for_row = ''
    for col in range(ncolumns):
        if col != 0:
            string_for_row += col_separation
        string_for_row += f'{{{col + 1}: >{{0}}}}'

    print(main_title)
    for (nrow, row) in enumerate(tabl_data):
        if nrow == 0:
            print((Format.underline + string_for_row + Format.end).format(char_in_column, *row))
        else:
            print(string_for_row.format(char_in_column, *row))


def standardize_dataset(input_spectra, input_magnitudes,
                        label_sfh, label_z,
                        method_input=1, method_label=None):

    if method_label is None:
        method_label = method_input

    # Inputs
    input_spectra_out = standardize_single_dataset(input_spectra, method_input)
    input_magnitudes_out = standardize_single_dataset(input_magnitudes, method_input)

    # Labels
    label_sfh_out = standardize_single_dataset(label_sfh, method_label)
    label_z_out = standardize_single_dataset(label_z, method_label)

    return input_spectra_out, input_magnitudes_out, label_sfh_out, label_z_out


def standardize_single_dataset(data, method):
    if method == 0:
        # Skip normalization
        output = data
    elif method == 1:
        # Normalize spectra for the middle_value = 1
        center_idx = int(len(data[0, :]) / 2 - 0.5)
        output = data / data[:, center_idx][:, None]
    elif method == 2:
        # Divide by the average
        output = data / data.mean(axis=0, keepdims=True)
    elif method == 3:
        # ToDo: Not implemented.
        # Stretch out and normalize to mean=1, std=1. (Maybe mean=0.5 or mean=0)
        Warning("Not implemented for method == 3")
        output = data
    else:
        raise ValueError("method used out of bounds.")
    return output
