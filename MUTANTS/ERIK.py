# -*- coding: utf-8 -*-
# filE ReadIng and Kleaning [Magneto]

import json
import os
import typing
import subprocess
import uuid
from shutil import copyfile
from datetime import datetime
import numpy as np
from astropy.io import fits

from RAVEN import standardize_single_dataset, convert_bytes


def loadfiles(input_path: str = "/Volumes/Elements/Outputs/Input_20211213T154548_HjCktf.fits",
              labels_path: str = "/Volumes/Elements/Outputs/Label_20211213T154548_HjCktf.fits",
              method_standardize_spectra=2,
              method_standardize_magnitudes=4,
              method_standardize_label_sfh=3,
              method_standardize_label_z=3,
              verbose=1) -> typing.Tuple[np.array, np.array, np.array, np.array, np.array, np.array]:
    """
    Load the dataset from file.

    :param input_path:
    :param labels_path:
    :param method_standardize_spectra:
    :param method_standardize_magnitudes:
    :param method_standardize_label_sfh:
    :param method_standardize_label_z:
    :param verbose:
    :return:
    """
    # ToDo: argparse these variables

    # Verify data is in 32 bits
    verify32bits(input_path)
    verify32bits(labels_path)

    print("[INFO] Loading inputs...")
    input_data, input_header = open_fits_file(input_path)
    print("[INFO] Loading label...")
    label_data, label_header = open_fits_file(labels_path)
    print("[INFO] Input Data shape:", input_data.shape)

    # Read the arrays from files
    input_lenght_spectra = input_header["nspectra"]
    label_agevec = label_header["nagevec"]
    spectra_lambda = input_data[1:input_lenght_spectra + 1, 0]
    agevec = label_data[1:label_agevec + 1, 0]
    input_spectra = input_data[1:input_lenght_spectra + 1, 1:]
    input_magnitudes = input_data[input_lenght_spectra + 1:, 1:]
    label_sfh = label_data[1:label_agevec + 1, 1:]
    label_z = label_data[label_agevec + 1:, 1:]
    # All these vectors are ordered the other way around. Transpose them
    input_spectra = input_spectra.transpose()
    input_magnitudes = input_magnitudes.transpose()
    label_sfh = label_sfh.transpose()
    label_z = label_z.transpose()

    input_spectra, input_magnitudes, label_sfh, label_z = \
        standardize_dataset(input_spectra, input_magnitudes, label_sfh, label_z,
                            method_standardize_spectra=method_standardize_spectra,
                            method_standardize_magnitudes=method_standardize_magnitudes,
                            method_standardize_label_sfh=method_standardize_label_sfh,
                            method_standardize_label_z=method_standardize_label_z)
    if verbose > 0:
        print(f"""
        Variable sizes:
            Input_spectra: {input_spectra.shape} - {convert_bytes(input_spectra.nbytes)}
            Input_magnitudes: {input_magnitudes.shape} - {convert_bytes(input_magnitudes.nbytes)}
            Label_sfh: {label_sfh.shape} - {convert_bytes(label_sfh.nbytes)}
            Label_z: {label_z.shape} - {convert_bytes(label_z.nbytes)}
        """)

    return input_spectra, input_magnitudes, label_sfh, label_z, spectra_lambda, agevec


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
                print(f"File ({os.path.basename(filepath)}) has BITPIX={hdul[0].header['bitpix']}.\n" +
                      f"Reducing to single precision (-32) ...")
            if verbose == 2:
                print(f"File ({filepath}) has BITPIX={hdul[0].header['bitpix']}.\n"
                      f"Reducing to single precision (-32) ...")
            # Reduce to -32
            hdul[0].data = hdul[0].data.astype(np.float32)
            hdul.flush()
            if verbose == 2:
                print(f"Correctly reduced {os.path.basename(filepath)} to -32.")


def clean_line(idx_line, line_str):
    line_str = line_str.strip()
    if line_str[:2] == "//" or len(line_str) == 0:
        return None

    if line_str.count("=") != 1:
        print(f"[ERROR WHILE READING CONFIG FILE] line {idx_line} does not have a '"
              f"'valid number of '='. Line skipped.")
        return None

    # Check if there is a comment behind a parameter
    if "//" in line_str:
        line_str = line_str.split("//")[0]

    split_line = line_str.split("=")
    split_line = [x.strip() for x in split_line]  # Clean spaces

    return tuple(split_line)


def getparametersfromid(filename, id_searched, verbose=0, returnfunction=False):
    """
    Returns the parameters that were used to generate a specific piece of information given an ID and the metadata file.
    It will verify if the the metadata has the "Combined" parameter, in which case it will seach for the subset id.
    :param verbose:
    :param filename:
    :param id_searched:
    :param returnfunction:
    :return:
    """

    # ToDo: Maybe verify UUIDs?
    # ToDo: Have an input that calls the R execution to generate the data again?
    # Open Metadata file
    with open(filename) as file_:
        data = json.load(file_)

    if verbose > 1:
        print(json.dumps(data, indent=4, sort_keys=True))
        pass

    # Verify if combined
    try:
        combined = data["Combined"]
    except KeyError:
        combined = False

    # If combined, look in which subset needs to be searched.
    # Additionally, look up which is the value that was added (and needs to be subtracted)
    id_subset = None
    if combined:
        for id_subset, last_id_this_subset in enumerate(data["Last_ID"]):
            if id_searched <= last_id_this_subset:
                break

        # Reduce id_searched according to the previous last_id, if id_subset > 0
        if id_subset > 0:
            id_searched -= data["Last_ID"][id_subset - 1]
        # Reduce data to subset for the search
        data = data[str(id_subset)]

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
        if returnfunction:
            final_dictionary["mfunction"] = mass_data_mfunc["func"]
            final_dictionary["zfunction"] = z_data_mfunc["func"]
        return final_dictionary


def read_config_file(filename, file_folder=None, reset_file=False, default_config_file="Data/default_config_file.txt"):
    # ToDo: Docstring
    if type(filename) is not str:
        raise KeyError("filename needs to be a string")
    if file_folder is None:
        file_folder = os.getcwd()

    # simplify full filename
    full_filename = os.path.join(file_folder, filename)

    # Check if file is there. If it isn't, generate a blank file with information templat
    if not os.path.isfile(full_filename) or reset_file:
        print(f"[INFO] There is no config file. A template will be created at {os.path.join(file_folder, filename)} .")
        copyfile(default_config_file, os.path.join(file_folder, filename))
        # ToDo: This path should be linked to the env_variable that reads where the data for the library is stored

    # Open file and read the configuration parameters.
    with open(full_filename, 'r') as f:
        lines = f.readlines()

    # Verify sintax and clean list
    cv_params = False
    cv_parameters = dict()
    parameters = dict()
    for idx, line in enumerate(lines):
        cleaned_line = clean_line(idx, line)
        if cleaned_line is not None:
            if cleaned_line[0] == "CVParams" and eval(cleaned_line[1]):
                cv_params = True
                continue
            if cv_params:
                tmp = eval(cleaned_line[1])
                if type(tmp) == tuple:
                    cv_parameters[cleaned_line[0]] = tmp[0]
                else:
                    cv_parameters[cleaned_line[0]] = tmp

            else:
                tmp = eval(cleaned_line[1])
                if type(tmp) == tuple:
                    parameters[cleaned_line[0]] = tmp[0]
                else:
                    parameters[cleaned_line[0]] = tmp

    return parameters, cv_parameters


def combine_datasets(file_list_sufixes: list, file_folder="", combined_output_sufix: str = "combined", overwrite=True,
                     whichrscript: str = "/usr/local/bin/Rscript"):
    """
    Combines n datasets into a single combined dataset, where n=len(file_list_sufixes).
    The files are located in file_folder (relative or absolute path). The outputs will be stored in the same folder.
    The three files will be called:
        -Input_[combined_output_sufix].fits
        -Label_[combined_output_sufix].fits
        -MetaD_[combined_output_sufix].fits

    :param file_list_sufixes:           Suffixes should not contain extensions
    :param file_folder:
    :param combined_output_sufix:
    :param overwrite:
    :param whichrscript:
    """

    # if the length of file_list_sufixes == ["all"], it will use all files.
    if file_list_sufixes == ["all"]:
        # Get all files
        file_list_sufixes = os.listdir(file_folder)
        # Filter only Inputs and keep only sufixes
        file_list_sufixes = [_[6:-5] for _ in file_list_sufixes if _[0:6] == "Input_"]
        # If combined is here, remove it
        file_list_sufixes = [_ for _ in file_list_sufixes if _ != "combined"]
        # Sort
        file_list_sufixes.sort()
        print(f"[INFO] All files in path are going to be used. Sufixes: {', '.join(file_list_sufixes)}")

    now = datetime.now()

    uuidi = str(uuid.uuid4())
    uuidl = str(uuid.uuid4())
    uuidm = str(uuid.uuid4())

    # Initialize variables
    metadata_combined = {"Combined": True, "Last_ID": []}
    last_id = None
    in_data = None
    in_header = None
    lab_data = None
    lab_header = None
    length_dataset = []
    for idx, file_prefix in enumerate(file_list_sufixes):
        print("[INFO] Reading prefix", file_prefix, "...")
        input_name = os.path.join(file_folder, "Input_" + file_prefix + ".fits")
        label_name = os.path.join(file_folder, "Label_" + file_prefix + ".fits")
        metadata_name = os.path.join(file_folder, "MetaD_" + file_prefix + ".rda")

        for file in [input_name, label_name, metadata_name]:
            # Verify file exists
            if not os.path.isfile(file):
                raise FileNotFoundError(f"File {file} not found.")

        # Read data and concatenate the information
        if idx == 0:
            # This is the first value. It will be used as a base
            # Input
            in_data, in_header = open_fits_file(input_name)
            in_header["filename"] = "Input_" + combined_output_sufix + ".fits"

            # Label
            lab_data, lab_header = open_fits_file(label_name)
            lab_header["filename"] = "Label_" + combined_output_sufix + ".fits"

            # UUIDs are no longer valid.
            for key, val in zip(["UUIDINP", "UUIDLAB", "UUIDMET"], [uuidi, uuidl, uuidm]):
                in_header[key] = val
                lab_header[key] = val

            # Add new keywords to header
            new_header_keywords = ["Datecomb",
                                   f"endrow{idx:02}",
                                   f"prefix{idx:02}"]
            new_header_values = [(now.strftime("%Y/%m/%d - %H:%M:%S"), "Date when combined"),
                                 (in_header["nrows"], f"Last value for the {idx}-th dataset"),
                                 (file_prefix, f"Prefix for the {idx}-th dataset")]
            for keyw, value in zip(new_header_keywords, new_header_values):
                in_header[keyw] = value
                lab_header[keyw] = value

            # Define last ID value
            last_id = in_data[0, -1]
            length_dataset.append(in_header["NROWS"])
        else:
            # Read new information
            tmp_in_data, tmp_in_header = open_fits_file(input_name)
            tmp_lab_data, tmp_lab_header = open_fits_file(label_name)

            # Update data
            tmp_in_data = tmp_in_data[:, 1:]
            tmp_lab_data = tmp_lab_data[:, 1:]
            # Update IDs
            tmp_in_data[0, :] = tmp_in_data[0, :] + last_id
            tmp_lab_data[0, :] = tmp_lab_data[0, :] + last_id
            last_id = tmp_in_data[0, -1]

            # Concatenate data with combined array
            in_data = np.concatenate([in_data, tmp_in_data], axis=1)
            lab_data = np.concatenate([lab_data, tmp_lab_data], axis=1)

            # Add new headers and keywords
            new_header_keywords = [f"endrow{idx:02}",
                                   f"prefix{idx:02}"]
            new_header_values = [(int(last_id), f"Last value for the {idx}-th dataset"),
                                 (file_prefix, f"Prefix for the {idx}-th dataset")]
            for keyw, value in zip(new_header_keywords, new_header_values):
                in_header[keyw] = value
                lab_header[keyw] = value
            length_dataset.append(tmp_in_header["NROWS"])

    # update nrows
    in_header["nrows"] = last_id
    lab_header["nrows"] = last_id

    # Combine Metadata
    # ToDo: read path to script from environment
    subprocess_str = " ".join(
        [whichrscript, "--vanilla", "/Users/enrique/Documents/GitHub/LOGAN-SFH/CombineMetadata.R",
         "uuidi=" + uuidi, "uuidl=" + uuidl, "uuidm=" + uuidm,
         "path=" + file_folder,
         "out=" + "MetaD_" + combined_output_sufix + ".rda",
         ] + file_list_sufixes)
    subprocess.call(subprocess_str, shell=True)

    # Save new files
    print(f"[INFO] Saving directory: {file_folder}")
    print(f"[INFO] Saving combined Input file to Input_{combined_output_sufix}.fits ...")
    fits.writeto(file_folder + "/Input_" + combined_output_sufix + ".fits", in_data, in_header, overwrite=overwrite)
    print(f"[INFO] Saving combined Label file to Label_{combined_output_sufix}.fits ...")
    fits.writeto(file_folder + "/Label_" + combined_output_sufix + ".fits", lab_data, lab_header, overwrite=overwrite)

    print("[INFO] Converting double precision to single precision...")
    verify32bits(file_folder + "/Label_" + combined_output_sufix + ".fits", verbose=0)
    verify32bits(file_folder + "/Input_" + combined_output_sufix + ".fits", verbose=0)

    print(f"[INFO] Number of lines: {' '.join([str(_) for _ in length_dataset])} -> Total: {int(in_header['nrows'])}")
    print(f"[INFO] Finished combining {len(file_list_sufixes)} files.")


def standardize_dataset(input_spectra, input_magnitudes,
                        label_sfh, label_z,
                        method_standardize_spectra=2,
                        method_standardize_magnitudes=4,
                        method_standardize_label_sfh=3,
                        method_standardize_label_z=3):
    # Inputs
    input_spectra_out, mean_spectra = standardize_single_dataset(input_spectra, method_standardize_spectra)
    input_magnitudes_out, _ = standardize_single_dataset(input_magnitudes, method_standardize_magnitudes, mean_spectra)

    # Labels
    label_sfh_out, _ = standardize_single_dataset(label_sfh, method_standardize_label_sfh)
    label_z_out, _ = standardize_single_dataset(label_z, method_standardize_label_z)

    return input_spectra_out, input_magnitudes_out, label_sfh_out, label_z_out


def open_fits_file(filename):
    if not os.path.isfile(filename):
        raise FileNotFoundError(f"File {filename} not found.")
    with fits.open(filename) as hdul:
        data = hdul[0].data
        header = hdul[0].header
    return data, header


def parse_argparse_config_file_default(args, default_config_file_path="Data/default_config_file.txt"):
    """
    Reads argparse and config_file. Keep based on priority: argparse > config_file > default.
    :param args:
    :param default_config_file_path:
    :return:
    """
    # Initialize
    def_parameters = None
    extra_data_from_default_config_file = False

    # Set argparse verbosity
    if args.verbose:
        args.verbose = 2
    elif args.quiet:
        args.verbose = 0
    else:
        args.verbose = 1
    del args.quiet

    # Read config_file if given. If it was not given, read the default parameters
    if args.config_file is not None:
        parameters, cv_parameters = read_config_file(args.config_file, reset_file=args.reset_config_file)
        def_parameters, _ = read_config_file(default_config_file_path, reset_file=False)
    else:
        parameters, cv_parameters = read_config_file(default_config_file_path, reset_file=False)

    # If a parameter was given through argparse, place it in the parameters list
    for keyword in list(vars(args).keys()):
        if vars(args)[keyword] is not None:
            parameters[keyword] = vars(args)[keyword]

    # If there is a file that was not given in parameter config, it will be read from the default config file.
    if args.config_file is not None:
        for keyword in def_parameters.keys():
            if (keyword not in parameters.keys()) and (keyword not in cv_parameters.keys()):
                parameters[keyword] = def_parameters[keyword]
                if not extra_data_from_default_config_file:
                    if parameters["verbose"] > 0:
                        print("[INFO] Missing parameters set to default. For more information, increase verbosity.")
                    extra_data_from_default_config_file = True
                if parameters["verbose"] > 1:
                    print(f"[INFO Parameters] Missing parameter {keyword} set to default value: {parameters[keyword]}")

    return parameters, cv_parameters


def prettyfy_json_file(filename, verbose=1, indent=4):
    """
    Given a path to a json file, prints the metadata, and saves the file again with indent=[4].
    :param filename:
    :param verbose:
    :param indent:
    """
    with open(filename, 'r+') as f:
        data = json.load(f)
        if verbose > 0:
            print(json.dumps(data, indent=indent))
        f.write(json.dumps(data, indent=indent))


if __name__ == "__main__":
    sufixes = ["20220427T192930_WViuSKxC", "20220428T130515_4QykcNJR", "20220428T154400_VxxFQPr2",
               "20220428T182828_dfn2FS3c", "20220503T140118_wYhTmpOD"]

    combine_datasets(["all"], file_folder="/Users/enrique/Documents/GitHub/LOGAN-SFH/KK", combined_output_sufix= "combined", overwrite = True)
    # combine_datasets(sufixes, file_folder="/Users/enrique/Documents/GitHub/LOGAN-SFH/KK", combined_output_sufix= "combined", overwrite = True)

