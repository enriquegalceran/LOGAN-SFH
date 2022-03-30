# -*- coding: utf-8 -*-
# Repository of Auxiliary scripts and Various spEcialized functioNs [Mystique]

import math
import numpy as np


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


def all_equal_ivo(lst):
    """
    Check if all elements are equal.
    Source: https://stackoverflow.com/questions/3844801/check-if-all-elements-in-a-list-are-identical
    :param lst:
    :return:
    """
    return not lst or lst.count(lst[0]) == len(lst)


def print_train_test_sizes(split1, main_title=""):
    """
    Outputs a table with the sizes of the different groups. Useful for debugging.
    :param split1:
    :param main_title:
    :return:
    """
    data = [
        [' ', 'train', 'test'],
        ['Spect', str(split1[0].shape), str(split1[1].shape)],
        ['Magnitudes', str(split1[2].shape), str(split1[3].shape)],
        ['SFH', str(split1[4].shape), str(split1[5].shape)],
        ['Metallicity', str(split1[6].shape), str(split1[7].shape)]
    ]
    print_table(data, main_title)


def print_build_model_parameters(spectr_arguments, magn_arguments, sfh_arguments, metal_arguments, main_title=""):
    """
    Print Parameters for every branch's parameters
    :param spectr_arguments:
    :param magn_arguments:
    :param sfh_arguments:
    :param metal_arguments:
    :param main_title:
    :return:
    """

    all_keywords = [*spectr_arguments.keys(), *magn_arguments.keys(),
                    *sfh_arguments.keys(), *metal_arguments.keys()]
    # Remove all repeated keywords
    all_keywords = list(set(all_keywords))
    data = [
        [' ', 'spectr', 'magn', 'sfh', 'metal']
    ]

    for keyword in all_keywords:
        tmp = [keyword]
        for branch in [spectr_arguments, magn_arguments, sfh_arguments, metal_arguments]:
            if keyword in branch.keys():
                tmp.append([branch[keyword]])
            else:
                tmp.append(["N/A"])
        data.append(tmp)

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


def standardize_single_dataset(data, method, input_mean_value=None):
    """
    Standardize the data depending on the method employed.
            :param data:        array-like of shape (n_samples, n_features)
            :param method:      int. Type of normalization that will be executed
                                0:  No normalization. output=input
                                1:  Normalize based on the center value of the data. Interesting for spectra.
                                    The center value taken is int(n_features / 2 - 0.5).
                                    This could generate some issues, as the center value of the wavelength will
                                    always be 1.
                                2:  Normalize for <data> = 1. This is the default value for the spectra.
                                    This will be the intended normalization for spectra.
                                3:  Stretch/shrink data to be between 0 and 1.
                                    output=(output-output.min)/output.max
                                    If a row is constant (i.e. constant metallicity), the return value will be the input
                                    value. This is intented for the labels (SFH and Z). If the metallicity is constant,
                                    it will already be set between [0, 1].
                                    This will be the intended normalization for the SFH and Metallicity.
                                4:  Normalization for magnitudes. Any other input data will not result in a valuable
                                    normalization.
                                    output=data-2.5*log10(input_mean_value)
                                    where input_mean_value is the mean value of the mean value of the spectra (saved
                                    from before).
                                    If input_mean_value is not given, it will be set to 1-> log(1)=0-> no normalization
                                5:  [Not Implemented yet] Normalizes the data so that the mean value equals X and
                                    standarddistribution equals Y.
                                    (Inspired in the gaussian normalization to mu=0, sigma=1)
                                    (Preemptive values for X and Y: 1, 1)
    :param input_mean_value:    Used for method==4. Mean value of the spectra, calculated previously.
    :return:
    """
    mean_value = data.mean(axis=1)
    if method == 0:
        # Skip normalization
        output = data

    elif method == 1:
        # Normalize spectra for the center wavelength = 1
        center_idx = int(len(data[0, :]) / 2 - 0.5)
        output = data / data[:, center_idx][:, None]

    elif method == 2:
        # Divide by the average
        output = data / data.mean(axis=1, keepdims=True)

    elif method == 3:
        # Stretch to: Minimum value->0 and maximum value -> 1
        data[1, :] = np.arange(data.shape[1])
        output = data.copy()
        output = output - output.min(axis=1, keepdims=True)

        # A problem arises if a row has a constant value. In that case, we would have operated 0/0->NaN.
        # To fix that, we will check which max value is zero and we will set the value to the original case.
        # For the rest, we will still operate normally.
        max_rows = output.max(axis=1)
        idx_constant = max_rows == 0.
        output[idx_constant] = data[idx_constant]
        output[np.invert(idx_constant)] = output[np.invert(idx_constant)] / max_rows[np.invert(idx_constant), None]

    elif method == 4:
        # For magnitudes
        if input_mean_value is None:
            Warning("No input has been given for mean value. will be set to 1 (data will be unchanged).")
            input_mean_value = np.ones((data.shape[0], ))
        # output: data - 2.5 * log10(input_mean_value) + 20
        # This "+20" is added systematically to every single value
        output = data + 2.5 * np.log10(input_mean_value)[:, None] + 20

    elif method == 5:
        # ToDo: Not implemented.
        # Stretch out and normalize to mean=1, std=1. (Maybe mean=0.5 or mean=0)
        Warning("Not yet implemented for method == 5. Used method == 1")
        output = standardize_single_dataset(data=data, method=1)

    else:
        raise ValueError("method used out of bounds [0-5].")

    return output, mean_value
