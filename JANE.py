# -*- coding: utf-8 -*-
# Just A Network Executor


# help("modules")
# help("modules tensorflow")
import typing
import pandas as pd

import numpy as np
from datetime import datetime
from astropy.io import fits
# import matplotlib.pyplot as plt
import os
import json
import random
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RandomizedSearchCV
from scikeras.wrappers import KerasRegressor

from TESSA import convert_bytes, print_train_test_sizes, standardize_dataset, open_fits_file
from XAVIER import Cerebro

print("[INFO] Finished Importing")


def loadfiles(input_path: str = "/Volumes/Elements/Outputs/Input_20211213T154548_HjCktf.fits",
              labels_path: str = "/Volumes/Elements/Outputs/Label_20211213T154548_HjCktf.fits",
              method_standardize_spectra=2,
              method_standardize_magnitudes=4,
              method__standardize_label=3) -> typing.Tuple[np.array, np.array, np.array, np.array, np.array, np.array]:
    """
    Load the dataset from file.

    :param input_path:
    :param labels_path:
    :param method_standardize_spectra:
    :param method_standardize_magnitudes:
    :param method__standardize_label:
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
    input_length_filters = input_header["nfilters"]
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
                            method_standardize_label=method__standardize_label)

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


def getparametersfromid(filename, id_searched, verbose=0):
    """
    Returns the parameters that were used to generate a specific piece of information given an ID and the metadata file.
    It will verify if the the metadata has the "Combined" parameter, in which case it will seach for the subset id.
    :param verbose:
    :param filename:
    :param id_searched:
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
        return final_dictionary


def read_config_file(filename, file_folder=None, reset_file=False):
    # Todo: Remove
    reset_file = True

    if type(filename) is not str:
        raise KeyError("filename needs to be a string")
    if file_folder is None:
        file_folder = os.getcwd()

    # simplify full filename
    full_filename = os.path.join(file_folder, filename)

    # Check if file is there. If it isn't, generate a blank file with information templat
    if not os.path.isfile(full_filename) or reset_file:
        print(f"[INFO] There is no config file. A template will be created at {os.path.join(file_folder, filename)} .")
        content_in_file = """
// Config File For JANE.py
// Instructions:
// Precede comments with two slashes ("//")
// Comments can be placed behind 


// General parameters
initial_learning_rate = 1e-3
train_size = 0.8
test_size = 0.20
data_path = "/Volumes/Elements/Outputs"
data_sufix = "20220209T142129_fJnn24"
output_model_path = "/Volumes/Elements/Outputs/"

// Prefixes for each branch of the neural network
spectra_prefix ="spect_"
magnitude_prefix = "magn_"
sfh_prefix = "sfh_"
metal_prefix = "metal_"


// Default parameters of the network.
// These can be changed, following the method of Cerebro.build_iterative_branch.
// (If a single value is present, each layer will have the same value. If the value is a list, each value will be used
// with the respective layer.)

// Spectra branch
spect_branch_type = "cnn"
spect_number_layers = 3
spect_neurons_first_layer = 128
spect_progression = 0.5
spect_filter_size = [30, 10, 3]
spect_stride = 1
spect_act = "relu"
spect_pool_size = 3
spect_dropout = 0.15
spect_explicit = False
spect_output = 1
spect_output_neurons = 128
spect_final_act = "relu"
spect_final_layer_name = "spectra_intermediate"
spect_kernel_initializer = "glorot_uniform"

// Magnitudes branch
magn_branch_type = "dense"
magn_number_layers = 2
magn_neurons_first_layer = 64
magn_progression = 0.5
magn_act = "relu"
magn_dropout = 0.25
magn_explicit = False,
magn_output = 1
magn_output_neurons = 32
magn_final_act = "relu"
magn_final_layer_name = "magnitude_intermediate"
magn_kernel_initializer = "glorot_uniform"

// Star Formation History Branch
sfh_branch_type = "dense"
sfh_layers = [512, 256, 256, 128]
sfh_act = "relu"
sfh_dropout = 0.15
sfh_explicit = False
sfh_output = 2
sfh_output_neurons = 17
sfh_final_act = "relu"
sfh_final_layer_name = "sfh_output"

// Metallicity Branch
metal_branch_type = "dense"
metal_layers = [512, 256, 256, 128]
metal_act = "relu"
metal_dropout = 0.15
metal_explicit = False
metal_output = 2
metal_output_neurons = 17
metal_final_act = "relu",
metal_final_layer_name = "metallicity_output"
metal_kernel_initializer = "glorot_uniform"

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Information regarding the construction of the CNN
// Information that will be passed to the Cross-Validation needs to be placed after "CVParams = True"
// (Uncomment next line)
// CVParams = True

// Place the parameters that will be iterated in a list as f
// spect_neurons_first_layer = 256

        """
        # ToDo: Finish this
        with open(full_filename, 'w+') as f:
            f.writelines(content_in_file)

    # Open file and read the configuration parameters.

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
        split_line = [x.strip() for x in split_line]    # Clean spaces

        return tuple(split_line)


    with open(full_filename, 'r') as f:
        lines = f.readlines()

    # Verify sintax and clean list
    parameters = dict()
    for idx, line in enumerate(lines):
        cleaned_line = clean_line(idx, line)
        if cleaned_line is not None:
            parameters[cleaned_line[0]] = eval(cleaned_line[1])


    return parameters




def combine_datasets(file_list_sufixes, file_folder="", combined_output_sufix="combined", overwrite=True):
    """
    Combines n datasets into a single combined dataset, where n=len(file_list_sufixes).
    The files are located in file_folder (relative or absolute path). The outputs will be stored in the same folder.
    The three files will be called:
        -Input_[combined_output_sufix].fits
        -Label_[combined_output_sufix].fits
        -MetaD_[combined_output_sufix].fits
    
    :param file_list_sufixes: 
    :param file_folder: 
    :param combined_output_sufix: 
    :param overwrite:
    """
    now = datetime.now()

    # Initialize variables
    metadata_combined = {"Combined": True, "Last_ID": []}
    last_id = None
    in_data = None
    in_header = None
    lab_data = None

    lab_header = None
    for idx, file_prefix in enumerate(file_list_sufixes):
        input_name = file_folder + "Input_" + file_prefix + ".fits"
        label_name = file_folder + "Label_" + file_prefix + ".fits"
        metadata_name = file_folder + "MetaD_" + file_prefix + ".json"

        for file in [input_name, label_name, metadata_name]:
            # Verify file exists
            if not os.path.isfile(file):
                raise FileNotFoundError(f"File {file} not found.")

        # Metadata
        with open(metadata_name) as f:
            metadata_combined[idx] = json.load(f)

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
            for key in ["UUIDINP", "UUIDLAB", "UUIDMET"]:
                del in_header[key]
                del lab_header[key]

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
            new_header_values = [(last_id, f"Last value for the {idx}-th dataset"),
                                 (file_prefix, f"Prefix for the {idx}-th dataset")]
            for keyw, value in zip(new_header_keywords, new_header_values):
                in_header[keyw] = value
                lab_header[keyw] = value

        # Add last_id to metadata
        metadata_combined["Last_ID"].append(int(last_id))

    # Save new files
    print(f"[INFO] Saving combined Input file to Input_{combined_output_sufix}.fits ...")
    fits.writeto(file_folder + "Input_" + combined_output_sufix + ".fits", in_data, in_header, overwrite=overwrite)
    print(f"[INFO] Saving combined Label file to Label_{combined_output_sufix}.fits ...")
    fits.writeto(file_folder + "Label_" + combined_output_sufix + ".fits", in_data, in_header, overwrite=overwrite)
    print(f"[INFO] Saving combined metadata file saving to MetaD_{combined_output_sufix}.fits ...")
    with open(file_folder + "MetaD_" + combined_output_sufix + ".json", 'w') as outfile:
        outfile.write(json.dumps(metadata_combined, indent=4))


def main(do_not_verify=True, **main_kwargs):
    # ToDo: Generate an argparser.

    # ToDo: Beautify the parameters
    # Parameters
    epochs = 50  # Number of epochs
    init_lr = 1e-3  # Initial learning rate
    bs = 32  # batches
    # Train, val and test sizes
    train_size = 0.80
    test_size = 0.20
    cv = 5
    random.seed(42)  # Set seed for testing purposes
    traintestrandomstate = 42  # Random state for train test split (default = None)
    traintestshuffle = True  # Shuffle data before splitting into train test (default = True)
    loss_function_used = "SMAPE"  # Define which lossfunction should be used (i.e. "SMAPE")
    # data_path = "/Volumes/Elements/Outputs/"
    data_path = ""
    data_sufix = "20220209T142129_fJnn24"
    output_model_path = "/Volumes/Elements/Outputs/"
    # path_output_plots = "/Volumes/Elements/Outputs/plot"
    output_plots_path = "plot"

    # Load Data
    print("[INFO] Loading data...")
    input_spectra, input_magnitudes, label_sfh, label_z, spectra_lambda, agevec = \
        loadfiles(input_path=data_path + "Input_" + data_sufix + ".fits",
                  labels_path=data_path + "Label_" + data_sufix + ".fits")
    # ToDo: Shuffle inputs and labels (together!!)

    print(f"""
    Variable sizes:
        Input_spectra: {input_spectra.shape} - {convert_bytes(input_spectra.nbytes)}
        Input_magnitudes: {input_magnitudes.shape} - {convert_bytes(input_magnitudes.nbytes)}
        Label_sfh: {label_sfh.shape} - {convert_bytes(label_sfh.nbytes)}
        Label_z: {label_z.shape} - {convert_bytes(label_z.nbytes)}
        """)

    # Split the data into training+validation and testing
    # ToDo: this section will be deprecated with the Cross-Validation
    assert train_size + test_size == 1, "The sum of train + test sizes has to add up to '1.0'."
    split_train_test = train_test_split(input_spectra, input_magnitudes, label_sfh, label_z,
                                        test_size=test_size,
                                        train_size=train_size,
                                        random_state=traintestrandomstate,
                                        shuffle=traintestshuffle)
    (trainSpect, testSpect,
     trainMag, testMag,
     trainLabSfh, testLabSfh,
     trainLabZ, testLabZ) = split_train_test

    # Concatenate inputs for single model
    trainData = np.concatenate([trainSpect, trainMag], axis=1)
    trainLabels = np.concatenate([trainLabSfh, trainLabZ], axis=1)
    testData = np.concatenate([testSpect, testMag], axis=1)
    testLabels = np.concatenate([testLabSfh, testLabZ], axis=1)

    # Verify the size of the
    print_train_test_sizes(split_train_test, main_title="Sizes of diferent sets (Train, Test)")
    # "split" tuple can be removed now
    del split_train_test

    ############################################################################
    # Build model
    custom_kwargs_model = {"spect_neurons_first_layer": 256}
    # model = Cerebro.build_model(epochs=epochs, loss_function_used=loss_function_used, init_lr=init_lr,
    #                             **custom_kwargs_model, **main_kwargs)
    # model.summary()
    # Cerebro.graph(model, "tstimage4.png")

    ############################################################################
    # Cross Validation
    # Define Regressor Model

    # Set parameters that will generate the grid
    # Here is where all the parameters will go (regarding 'spect_', ...)
    param_grid = dict(
        epochs=[50, 75, 100],
        batch_size=[50, 100, 200],
        magn_neurons_first_layer=[32, 64, 128],
        magn_number_layers=[2, 3],
        spect_number_layers=[4],
        spect_neurons_first_layer=[128, 256],
        spect_filter_size=[[30, 20, 10, 5], [30, 30, 10, 5], [50, 25, 10, 5], [30, 15, 10, 5], 15],
    )
    print("param_grid", param_grid)
    number_of_combinations = 1
    for key, value in param_grid.items():
        print(f"{key}: {value} - {len(value)}")
        number_of_combinations *= len(value)
    print(f"[INFO] Number of possible combinations: {number_of_combinations}")

    estimator = KerasRegressor(model=Cerebro.build_model, verbose=1, **param_grid)
    print("estimator", estimator.get_params().keys())

    grid = RandomizedSearchCV(estimator=estimator, param_distributions=param_grid,
                              n_iter=10,
                              cv=cv, random_state=traintestrandomstate, verbose=1)
    # grid = GridSearchCV(estimator=estimator, param_grid=param_grid, cv=cv, verbose=5)
    # print("grid", grid.get_params().keys())

    # Train model
    if not do_not_verify:
        continue_with_training = input("Continue training? [Y]/N ")
        if continue_with_training.lower() in ["n", "no", "stop"]:
            raise KeyboardInterrupt('User stopped the execution.')
    print("[INFO] Start training...")
    grid_result = grid.fit(trainData, trainLabels)

    # grid_result = grid.fit(x={"spectra_input": trainSpect, "magnitude_input": trainMag},
    #                        y={"sfh_output": trainLabSfh, "metallicity_output": trainLabZ})

    # print results
    print("\n\n\n")
    print(grid_result)

    print(f'Best Accuracy for {grid_result.best_score_:.4} using {grid_result.best_params_}')
    means = grid_result.cv_results_['mean_test_score']
    stds = grid_result.cv_results_['std_test_score']
    params = grid_result.cv_results_['params']
    for mean, stdev, param in zip(means, stds, params):
        print(f'mean={mean:.4}, std={stdev:.4} using {param}')

    accuracy = grid.score(testData, testLabels)
    pd.set_option('display.max_columns', None)
    print(f"\n\n\nThe test accuracy score of the best model is "
          f"{accuracy:.2f}")
    from pprint import pprint

    print("The best parameters are:")
    pprint(grid.best_params_)
    # get the parameter names
    column_results = [
        f"param_{name}" for name in param_grid.keys()]
    column_results += [
        "mean_test_score", "std_test_score", "rank_test_score"]

    cv_results = pd.DataFrame(grid.cv_results_)
    cv_results = cv_results[column_results].sort_values(
        "mean_test_score", ascending=False)
    print(cv_results)

    cv_results = cv_results.set_index("rank_test_score")
    print(cv_results["mean_test_score"][1] - cv_results["mean_test_score"][2])

    # train_history = model.fit(x={"spectra_input": trainSpect, "magnitude_input": trainMag},
    #                           y={"sfh_output": trainLabSfh, "metallicity_output": trainLabZ},
    #                           validation_data=({"spectra_input": testSpect, "magnitude_input": testMag},
    #                                            {"sfh_output": testLabSfh, "metallicity_output": testLabZ}),
    #                           epochs=epochs,
    #                           batch_size=bs,
    #                           verbose=1)
    #
    # # save the model to disk
    # print("[INFO] serializing network...")
    # path_model = path_output_model_path + "/test_model_" + data_sufix + "_" + datetime.now().strftime("%d%m%YT%H%M%S") + ".h5"
    # print("[INFO] Model stored in", path_model)
    # model.save(path_model, save_format="h5")

    # ############
    # # ToDo: Consolidate into a single function
    # # plot the total loss, category loss, and color loss
    # loss_names = ["loss", "sfh_output_loss", "metallicity_output_loss"]
    # plt.style.use("ggplot")
    # (fig, ax) = plt.subplots(3, 1, figsize=(13, 13))
    # # loop over the loss names
    # for (i, l) in enumerate(loss_names):
    #     # plot the loss for both the training and validation data
    #     title = "Loss for {}".format(l) if l != "loss" else "Total loss"
    #     ax[i].set_title(title)
    #     ax[i].set_xlabel("Epoch #")
    #     ax[i].set_ylabel("Loss")
    #     ax[i].plot(np.arange(0, epochs), train_history.history[l], label=l)
    #     ax[i].plot(np.arange(0, epochs), train_history.history["val_" + l], label="val_" + l)
    #     ax[i].legend()
    # # save the losses figure
    # plt.tight_layout()
    # plt.savefig("{}_losses.png".format(path_output_plots))
    # print("[INFO] Loss image stored in {}_losses.png".format(path_output_plots))
    # plt.close()
    #
    # # create a new figure for the accuracies
    # accuracy_names = ["sfh_output_accuracy", "metallicity_output_accuracy"]
    # plt.style.use("ggplot")
    # (fig, ax) = plt.subplots(2, 1, figsize=(8, 8))
    # # loop over the accuracy names
    # for (i, l) in enumerate(accuracy_names):
    #     # plot the loss for both the training and validation data
    #     ax[i].set_title("Accuracy for {}".format(l))
    #     ax[i].set_xlabel("Epoch #")
    #     ax[i].set_ylabel("Accuracy")
    #     ax[i].plot(np.arange(0, epochs), train_history.history[l], label=l)
    #     ax[i].plot(np.arange(0, epochs), train_history.history["val_" + l], label="val_" + l)
    #     ax[i].legend()
    # # save the accuracies figure
    # plt.tight_layout()
    # plt.savefig("{}_accs.png".format(path_output_plots))
    # print("[INFO] Acc image stored in {}_accs.png".format(path_output_plots))
    # plt.close()

    print("[INFO] Finished!")


if __name__ == "__main__":
    read_config_file("config.txt")
    # main()
