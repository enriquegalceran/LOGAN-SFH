# -*- coding: utf-8 -*-
# Just A Network Executor
import typing

import numpy as np
from datetime import datetime
from astropy.io import fits
import matplotlib.pyplot as plt
import os
import json
import random

from TESSA import convert_bytes, print_train_test_val_sizes, standardize_dataset
from XAVIER import Cerebro

# For execution
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import train_test_split


def loadfiles(input_path: str = "/Volumes/Elements/Outputs/Input_20211213T154548_HjCktf.fits",
              labels_path: str = "/Volumes/Elements/Outputs/Label_20211213T154548_HjCktf.fits",
              size_inputs=None,
              size_magnitudes=None,
              size_labels=None,
              method_input=1,
              method_label=None) -> typing.Tuple[np.array, np.array, np.array, np.array, np.array, np.array]:
    """
    Load the dataset from file.

    :param input_path:
    :param labels_path:
    :param size_inputs:
    :param size_magnitudes:
    :param size_labels:
    :param method_input:
    :param method_label:
    :return:
    """
    # ToDo: argparse these variables

    # Verify data is in 32 bits
    verify32bits(input_path)
    verify32bits(labels_path)

    print("[INFO] Loading inputs...")
    with fits.open(input_path) as hdul:
        input_data = hdul[0].data
        input_header = hdul[0].header

    print("[INFO] Loading label...")
    with fits.open(labels_path) as hdul:
        label_data = hdul[0].data
        label_header = hdul[0].header

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
                            method_input=method_input, method_label=method_label)

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
    # ToDo: Generate an argparser.

    # ToDo: Beautify the parameters
    # Parameters
    epochs = 50  # Number of epochs
    init_lr = 1e-3  # Initial learning rate
    bs = 32  # batches
    # Train, val and test sizes
    train_size = 0.70
    val_size = 0.10
    test_size = 0.20
    random.seed(42)  # Set seed for testing purposes
    traintestrandomstate = 42  # Random state for train test split (default = None)
    traintestshuffle = True  # Shuffle data before splitting into train test (default = True)
    loss_function_used = "SMAPE"  # Define which lossfunction should be used ("crossentropy"/"SMAPE")
    # data_path = "/Volumes/Elements/Outputs/"
    data_path = ""
    data_sufix = "20220119T154253_Hr3TXx"
    path_output_model_path = "/Volumes/Elements/Outputs/"
    # path_output_plots = "/Volumes/Elements/Outputs/plot"
    path_output_plots = "plot"

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
    assert train_size + val_size + test_size == 1, "The sum of the three train/val/test sizes has to add up to '1.0'."
    train_val_size = train_size + val_size
    split_trainval_test = train_test_split(input_spectra, input_magnitudes, label_sfh, label_z,
                                           test_size=test_size,
                                           train_size=train_val_size,
                                           random_state=traintestrandomstate,
                                           shuffle=traintestshuffle)
    (trainvalSpect, testSpect,
     trainvalMag, testMag,
     trainvalLabSfh, testLabSfh,
     trainvalLabZ, testLabZ) = split_trainval_test

    # Split the training+validation into training and validation
    split_train_val = train_test_split(trainvalSpect, trainvalMag, trainvalLabSfh, trainvalLabZ,
                                       test_size=val_size / train_val_size,
                                       train_size=train_size / train_val_size,
                                       random_state=traintestrandomstate,
                                       shuffle=traintestshuffle)
    (trainSpect, valSpect,
     trainMag, valMag,
     trainLabSfh, valLabSfh,
     trainLabZ, valLabZ) = split_train_val

    # Verify the size of the
    print_train_test_val_sizes(split_trainval_test, split_train_val,
                               main_title="Sizes of diferent sets (Train, Val, Test)")
    # "split" tuples can be removed now
    del split_trainval_test
    del split_train_val

    # Build model
    print("[INFO] Building model...")
    model = Cerebro.build_model(spectra_data_shape=3761, magnitudes_data_shape=5,
                                number_neurons_spec=256, number_neurons_magn=32,
                                number_output_sfh=trainLabSfh.shape[1], number_output_metal=trainLabZ.shape[1],
                                explicit=False)
    model.summary()
    # Cerebro.graph(model, "tstimage.png")

    # Define loss. Can be different functions. Right now, only crossentropy and smape.
    if loss_function_used == "crossentropy":
        # Standard crossentropy
        losses = {
            "sfh_output": "categorical_crossentropy",
            "metallicity_output": "categorical_crossentropy"
        }
    elif loss_function_used == "SMAPE":
        # SMAPE
        # Note: indicate the function, do not call it. (i.e.: 'Cerebro.smape_loss'; NOT 'Crebro.smape_loss()')
        losses = {
            "sfh_output": Cerebro.smape_loss,
            "metallicity_output": Cerebro.smape_loss
        }
    else:
        raise ValueError("loss_function_used not defined/not one of the accepted values!")

    # Loss weight for the two different branches
    loss_weights = {
        "sfh_output": 1.0,
        "metallicity_output": 0.8
    }

    # Initialize optimizer and compile the model
    print("[INFO] Compiling model...")
    opt = Adam(lr=init_lr, decay=init_lr / epochs)
    model.compile(optimizer=opt, loss=losses, loss_weights=loss_weights, metrics=["accuracy"])

    # Train model
    continue_with_training = input("Continue training? [Y]/N")
    if continue_with_training in ["N", "no", "stop"]:
        exit(1)

    print("[INFO] Start training...")
    train_history = model.fit(x={"spectra_input": trainSpect, "magnitude_input": trainMag},
                              y={"sfh_output": trainLabSfh, "metallicity_output": trainLabZ},
                              validation_data=({"spectra_input": testSpect, "magnitude_input": testMag},
                                               {"sfh_output": testLabSfh, "metallicity_output": testLabZ}),
                              epochs=epochs,
                              batch_size=bs,
                              verbose=1)

    # save the model to disk
    print("[INFO] serializing network...")
    path_model = path_output_model_path + "/test_model_" + data_sufix + "_" + datetime.now().strftime("%d%m%YT%H%M%S") + ".h5"
    print("[INFO] Model stored in", path_model)
    model.save(path_model, save_format="h5")

    ############
    # ToDo: Consolidate into a single function
    # plot the total loss, category loss, and color loss
    loss_names = ["loss", "sfh_output_loss", "metallicity_output_loss"]
    plt.style.use("ggplot")
    (fig, ax) = plt.subplots(3, 1, figsize=(13, 13))
    # loop over the loss names
    for (i, l) in enumerate(loss_names):
        # plot the loss for both the training and validation data
        title = "Loss for {}".format(l) if l != "loss" else "Total loss"
        ax[i].set_title(title)
        ax[i].set_xlabel("Epoch #")
        ax[i].set_ylabel("Loss")
        ax[i].plot(np.arange(0, epochs), train_history.history[l], label=l)
        ax[i].plot(np.arange(0, epochs), train_history.history["val_" + l], label="val_" + l)
        ax[i].legend()
    # save the losses figure
    plt.tight_layout()
    plt.savefig("{}_losses.png".format(path_output_plots))
    print("[INFO] Loss image stored in {}_losses.png".format(path_output_plots))
    plt.close()

    # create a new figure for the accuracies
    accuracy_names = ["sfh_output_accuracy", "metallicity_output_accuracy"]
    plt.style.use("ggplot")
    (fig, ax) = plt.subplots(2, 1, figsize=(8, 8))
    # loop over the accuracy names
    for (i, l) in enumerate(accuracy_names):
        # plot the loss for both the training and validation data
        ax[i].set_title("Accuracy for {}".format(l))
        ax[i].set_xlabel("Epoch #")
        ax[i].set_ylabel("Accuracy")
        ax[i].plot(np.arange(0, epochs), train_history.history[l], label=l)
        ax[i].plot(np.arange(0, epochs), train_history.history["val_" + l], label="val_" + l)
        ax[i].legend()
    # save the accuracies figure
    plt.tight_layout()
    plt.savefig("{}_accs.png".format(path_output_plots))
    print("[INFO] Acc image stored in {}_accs.png".format(path_output_plots))
    plt.close()

    print("[INFO] Finished!")


if __name__ == "__main__":
    main()
