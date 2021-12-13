# Just A Network Executor

import numpy as np
import pandas as pd
import argparse
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import json


def loadfiles(input_path="/Volumes/Elements/Outputs/Input_20211213T115406_akvGb3.fits",
              labels_path="/Volumes/Elements/Outputs/Label_20211213T115406_akvGb3.fits"):
    # ToDo: argparse these variables

    # Verify data is in 32 bits
    verify32bits(input_path)
    verify32bits(labels_path)

    print("Loading inputs...")
    with fits.open(input_path) as hdul:
        input_data = hdul[0].data
        input_header = hdul[0].header

    print("Loading label...")
    with fits.open(labels_path) as hdul:
        label_data = hdul[0].data
        label_header = hdul[0].header

    print(input_data.shape)


def verify32bits(filepath, verbose=1):
    """
    Verifies if the file is in -32 or -64 format. If double (-64), reduces to single (-32).
    :param verbose:
    :param filepath:
    :return:
    """
    with fits.open(filepath, mode="update") as hdul:
        if hdul[0].header["bitpix"] == -64:
            if verbose == 1:
                print(f"File ({os.path.basename(filepath)}) has BITPIX={hdul[0].header['bitpix']}. \
                        Reducing to single precision (-32) ...")
            if verbose == 2:
                print(f"File ({filepath}) has BITPIX={hdul[0].header['bitpix']}. \
                        Reducing to single precision (-32) ...")
            # Reduce to -32
            hdul[0].data = hdul[0].data.astype(np.float32)
            hdul.flush()
            if verbose == 2:
                print(f"Correctly reduced {os.path.basename(filepath)} to -32.")


def getparametersfromid(filename, id_searched, verbose=0):
    """
    Returns the parameters that were used to generate a specific piece of information given an ID and the metadata file.
    :param verbose:
    :param filename:
    :param id_searched:
    :return:
    """

    # ToDo: Maybe verify UUIDs?
    # ToDo: Have an input that calls the R execution to generate the data again?
    # Open Metadata file
    with open(filename) as f:
        data = json.load(f)

    if verbose > 1:
        print(json.dumps(data, indent=4, sort_keys=True))
        pass

    # Read parameters that will be used
    random_samples = data["randomSamples"][0]
    order_parameters = data["orderParameters"]
    massfunc_names = list(order_parameters.keys())
    accumulated_combinations = 0

    # Iterate over the different massfunc Names
    for mfunc in massfunc_names:
        # Mass data for mfunc
        mass_data_mfunc = data["massParams"][mfunc]
        # Name of parameters
        mass_keys_mfunc = [x if x in list(mass_data_mfunc.keys()) else None for x in order_parameters[mfunc]["mass"]]
        # Possible values of parameters
        mass_parameters = [mass_data_mfunc[x] for x in mass_keys_mfunc]
        # Number of possible values for each parameters
        number_values_mass_parameters = [len(x) for x in mass_parameters]

        # Obtain same values for Z
        z_data_mfunc = data["ZParams"]
        z_keys_for_mfunc = [x if x in list(z_data_mfunc.keys()) else None for x in order_parameters[mfunc]["Z"]]
        z_parameters = [z_data_mfunc[x] for x in z_keys_for_mfunc]
        number_values_z_parameters = [len(x) for x in z_parameters]

        # Once all the data is recollected, number of cases are calculated
        # All the parameter names
        all_parameters = mass_keys_mfunc + z_keys_for_mfunc
        # Values of the parameters
        values_all_parameters = mass_parameters + z_parameters
        # How many parameters are there (+ randomSample)
        nparam = len(all_parameters) + 1
        # How many cases are there for each parameter
        number_all_parameters = number_values_mass_parameters + number_values_z_parameters + [random_samples]

        # Calculate how many iterations there are for every case
        number_combinations = [0] * nparam
        number_combinations[-1] = random_samples + 1
        for i in reversed(range(nparam - 1)):
            number_combinations[i] = number_combinations[i + 1] * number_all_parameters[i]

        # Verify if ID is bigger than all possible combinations for this massfunc
        # If true, skip current massfunc and try with the next. Increase accumulated_combinations
        if id_searched > accumulated_combinations + number_combinations[0]:
            accumulated_combinations += number_combinations[0]
            continue

        # If smaller, it will stay with this massfunc
        current_id = id_searched - accumulated_combinations - 1
        idx_param = [0] * nparam
        for idx in range(nparam - 1):
            # Calculate from biggest to smallest the index of the parameter that was used.
            idx_param[idx] = int(current_id / number_combinations[idx + 1])
            current_id -= idx_param[idx] * number_combinations[idx + 1]
        # Add randomSample at the end
        idx_param[-1] = current_id

        # Generate the final dictionary that will be returned
        final_dictionary = {"massfunction": mfunc}
        for f in range(nparam - 1):
            final_dictionary[all_parameters[f]] = values_all_parameters[f][idx_param[f]]
        if idx_param[-1] == 0:
            final_dictionary["randomSample"] = False
        else:
            final_dictionary["randomSample"] = True

        if verbose >= 1:
            print(final_dictionary)
        return final_dictionary


def main():
    # Get the parameter with which the data was generated from ID and metadatafile
    # for id_ in range(1, 73):
    #     getparametersfromid("MetadataOutput.json", id_, verbose=1)

    loadfiles()


if __name__ == "__main__":
    main()
