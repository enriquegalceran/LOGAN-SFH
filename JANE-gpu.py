# -*- coding: utf-8 -*-
# Just A Network Executor [Phoenix]


import argparse
import random
import sys
from pprint import pprint
print(sys.path)
# ToDo: This is a hack, but it works...
# Add MUTANTS to path
sys.path.append("/home/egalcera/github/LOGAN-SFH/MUTANTS")


# import joblib
import tensorflow as tf
import keras
import matplotlib.pyplot as plt
# help("modules")
# help("modules tensorflow")
import pandas as pd
from keras.callbacks import EarlyStopping, ModelCheckpoint
#from scikeras.wrappers import KerasRegressor
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import train_test_split

from MUTANTS.ERIK import *
from MUTANTS.RAVEN import convert_bytes, print_train_test_sizes
from MUTANTS.XAVIER import Cerebro

print("[INFO] Finished Importing")


# Print GPUs available
print("Num GPUs Available: ", len(tf.config.list_physical_devices('GPU')))
print("GPUs:", tf.config.list_physical_devices('GPU'))

# Log where the jobs are being executed
# This line should be used with care, as it fills the log
#tf.debugging.set_log_device_placement(True)


def main(**main_kwargs):

    # Argument Parser (argparse)
    parser = argparse.ArgumentParser(description="Generate and train a Neural Network")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true",
                       help="Increase Verbosity.")
    group.add_argument("-q", "--quiet", action="store_true",
                       help="Decrease Verbosity.")
    parser.add_argument("--no-confirm", action="store_true", default=None,
                        help="Use to skip confirmation before training.")
    parser.add_argument("-cf", "--config-file", default=None, type=str,
                        help="Path to config file")
    parser.add_argument("-rcf", "--reset-config-file", action="store_true", default=None,
                        help="Resets config file back to the default value.")
    parser.add_argument("-dp", "--data-path", default=None, type=str,
                        help="Path to Data.")
    parser.add_argument("-ds", "--data-sufix", default=None, type=str,
                        help="Data sufix.")
    parser.add_argument("-cv", default=None, type=int,
                        help="Number of folds for cross-validation")
    parser.add_argument("--train-size", default=None, type=float,
                        help="Train size. Needs to be in the interval (0,1)")
    parser.add_argument("-e", "--epochs", default=None, type=int,
                        help="Set number of epochs.")
    parser.add_argument("--test-size", default=None, type=float,
                        help="Test size. Needs to be in the interval (0,1)")
    parser.add_argument("--combine-datasets", default=None, type=str, nargs="+",
                        help="Combines multiple datasets into a single one. Give multiple sufixes separated by spaces.")
    parser.add_argument("--combine-outname", default="combined", type=str,
                        help="Name for the combined datasets. Requires the use of combine-datasets.")
    parser.add_argument("--combine-path", default=None, type=str,
                        help="Path to datasets to be combined. Requires the use of combine-datasets.")
    args = parser.parse_args()
    if args.reset_config_file and args.config_file is None:
        parser.error("--reset_config_file requires --config_file.")

    # Verify if datasets should be combined
    if args.combine_path is not None and args.combine_datasets is None:
        parser.error("--combine-path requires --combine-datasets.")
    if args.combine_datasets is not None:
        if len(args.combine_datasets) < 2:
            parser.error("--combine-datasets needs to have at least 2 elements.")
        else:
            if args.combine_path is None:
                args.combine_path = os.getcwd()
            if not os.path.isabs(args.combine_path):
                args.combine_path = os.path.join(os.getcwd(), args.combine_path)
            print("[INFO] Combining datasets:")
            print("Path:", args.combine_path)
            print("Datasets:", args.combine_datasets)
            print("Outputname:", args.combine_outname)
            combine_datasets(args.combine_datasets, file_folder=args.combine_path,
                             combined_output_sufix=args.combine_outname, overwrite=True)
            print("[INFO] Skipping rest of JANE.")
            return 0

    # ToDo: Remove next line
    if "config_file_path" in main_kwargs.keys():
        args.config_file = main_kwargs["config_file_path"]
    if "no_confirm" in main_kwargs.keys():
        args.no_confirm = main_kwargs["no_confirm"]

    # Parameters
    parameters, cv_parameters = parse_argparse_config_file_default(args)
    train_size = parameters["train_size"]
    test_size = parameters["test_size"]
    cv = parameters["cv"]
    random.seed(parameters["random_seed"])  # Set seed for testing purposes
    traintestrandomstate = parameters["random_seed"]  # Random state for train test split (default = None)
    traintestshuffle = parameters["traintestshuffle"]  # Shuffle data before splitting into train test (default = True)
    # loss_function_used = parameters["loss_function_used"]  # Define which lossfunction should be used (i.e. "SMAPE")
    # data_path = "/Volumes/Elements/Outputs/"
    # ToDo: Fix this in the final version (i.e. no more debug)
    data_path = parameters["data_path"]
    # data_path = ""
    #data_path = "/Volumes/Elements/Outputs/"
    # parameters["output_model_path"] = ""
    # parameters["output_model_path"] = "/Volumes/Elements/Outputs/"
    data_sufix = parameters["data_sufix"]
    # data_sufix = "combined"
    output_model_path = parameters["output_model_path"]
    # path_output_plots = "/Volumes/Elements/Outputs/plot"
    # output_plots_path = "plot"

    # Load Data
    print("[INFO] Loading data...")
    # Verify if there are indications for standardizing data:
    arguments_standardize = {}
    for key in parameters.keys():
        if "method_standardize_" in key:
            arguments_standardize[key] = parameters[key]
    print("Arguments Standardization")
    print(arguments_standardize)
    input_spectra, input_magnitudes, label_sfh, label_z, spectra_lambda, agevec, ageweights = \
        loadfiles(input_path=data_path + "Input_" + data_sufix + ".fits",
                  labels_path=data_path + "Label_" + data_sufix + ".fits",
                  **arguments_standardize)

    # Split the data into training+validation and testing
    assert 0 < train_size < 1, "train_size needs to be in the interval (0,1)."
    assert 0 < test_size < 1, "test_size needs to be in the interval (0,1)."
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
    if parameters["input_mode"] == "single":
        train_data = np.concatenate([trainSpect, trainMag], axis=1)
        test_data = np.concatenate([testSpect, testMag], axis=1)
    else:
        train_data = (trainSpect, trainMag)
        test_data = (testSpect, testMag)
    if parameters["output_mode"] == "single":
        train_labels = np.concatenate([trainLabSfh, trainLabZ], axis=1)
        test_labels = np.concatenate([testLabSfh, testLabZ], axis=1)
    else:
        train_labels = (trainLabSfh, trainLabZ)
        test_labels = (testLabSfh, testLabZ)

    # Verify the size of the
    print_train_test_sizes(split_train_test, main_title="Sizes of diferent sets (Train, Test)")
    # "split" tuple can be removed now
    del split_train_test

    ############################################################################
    # Build model
    model = Cerebro.build_model(**main_kwargs, **parameters)
    model.summary()
    # Cerebro.graph(model, "tstimage.png")

    # Callbacks
    best_model_path = os.path.join(parameters["output_model_path"], parameters["output_name_best_estimator"])
    es = EarlyStopping(monitor="val_loss", mode="min", verbose=1, patience=parameters["patience"])
    mc = ModelCheckpoint(best_model_path, monitor="val_loss", mode="min", save_best_only=True, verbose=1)
    print(f"[INFO] Best Model saved in: {best_model_path}")
    print(f"[INFO] Number of epochs: {parameters['epochs']}")

    # Train
    if not parameters["no_confirm"]:
        continue_with_training = input("Continue training? [Y]/N ")
        if continue_with_training.lower() in ["n", "no", "stop"]:
            raise KeyboardInterrupt('User stopped the execution.')
    print("[INFO] Start training...")
    history = model.fit(train_data, train_labels, validation_data=(test_data, test_labels),
                        batch_size=parameters["batch_size"],
                        epochs=parameters["epochs"],
                        callbacks=[es, mc])

    # evaluate the model
    # print("[Evaluate] Evaluating final model")
    # _, train_acc = model.evaluate(train_data, train_labels, verbose=1)
    # _, test_acc = model.evaluate(test_data, test_labels, verbose=1)
    # print('Train: %.3f, Test: %.3f' % (train_acc, test_acc))

    epochs_draw = len(history.epoch)
    if parameters["output_mode"] == "single":
        # plot training history
        plt.style.use("ggplot")
        (fig, ax) = plt.subplots(2, 1, figsize=(13, 13))

        # Loss
        loss_names = ["loss", "val_loss"]
        ax[0].set_title("Loss")
        ax[0].set_xlabel("Epoch #")
        ax[0].set_ylabel("Loss")
        for lo in loss_names:
            ax[0].plot(np.arange(0, epochs_draw), history.history[lo], label=lo)
        ax[0].legend()

        # Accuracy
        acc_names = ["accuracy", "val_accuracy"]
        ax[1].set_title("Accuracy")
        ax[1].set_xlabel("Epoch #")
        ax[1].set_ylabel("Loss")
        for ac in acc_names:
            ax[1].plot(np.arange(0, epochs_draw), history.history[ac], label=ac)
        ax[1].legend()

        # plt.show()
        # save the losses figure
        # plt.tight_layout()
        # plt.savefig("{}_losses.png".format(path_output_plots))
        # print("[INFO] Loss image stored in {}_losses.png".format(path_output_plots))
        # plt.close()

    elif parameters["output_mode"] == "double":
        # # plot the total loss, category loss, and color loss
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
            ax[i].semilogy(np.arange(0, epochs_draw), history.history[l], label=l)
            ax[i].semilogy(np.arange(0, epochs_draw), history.history["val_" + l], label="val_" + l)
            ax[i].legend()
        # # save the losses figure
        # plt.tight_layout()
        # plt.savefig("{}_losses.png".format(path_output_plots))
        # print("[INFO] Loss image stored in {}_losses.png".format(path_output_plots))
        # plt.close()
        plt.savefig("imagen1.png")
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
            ax[i].plot(np.arange(0, epochs_draw), history.history[l], label=l)
            ax[i].plot(np.arange(0, epochs_draw), history.history["val_" + l], label="val_" + l)
            ax[i].legend()
        # # save the accuracies figure
        # plt.tight_layout()
        # plt.savefig("{}_accs.png".format(path_output_plots))
        # print("[INFO] Acc image stored in {}_accs.png".format(path_output_plots))
        # plt.close()
        plt.savefig("imagen2.png")


    # load the saved model
    # saved_model = keras.models.load_model(best_model_path, custom_objects={"smape_loss": Cerebro.smape_loss})
    # evaluate the model
    print("[Evaluate] Evaluating best model")
    print("skipping")
    # _, train_acc = saved_model.evaluate(train_data, train_labels, verbose=1)
    # _, test_acc = saved_model.evaluate(test_data, test_labels, verbose=1)
    # print('Train: %.3f, Test: %.3f' % (train_acc, test_acc))









    plt.savefig("imagen.png")
    print("[INFO] Finished!")


if __name__ == "__main__":
    # print(getparametersfromid("/Volumes/Elements/Outputs/MetaD_combined.json", 34509, verbose=0))
    # combine_datasets(["20220222T112103_drDhe7", "20220222T134553_qbkxpF", "20220222T213001_nHS1Bf"],
    #                  "/Volumes/Elements/Outputs/",
    #                  )
    main()
