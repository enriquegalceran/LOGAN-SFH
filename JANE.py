# -*- coding: utf-8 -*-
# Just A Network Executor [Phoenix]


# help("modules")
# help("modules tensorflow")
import pandas as pd
import numpy as np
import argparse
# import matplotlib.pyplot as plt
import random
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RandomizedSearchCV
from scikeras.wrappers import KerasRegressor
# from keras.callbacks import EarlyStopping, ModelCheckpoint

from ERIK import loadfiles, parse_argparse_config_file_default
from RAVEN import convert_bytes, print_train_test_sizes
from XAVIER import Cerebro

print("[INFO] Finished Importing")


def main(do_not_verify=True, **main_kwargs):

    # Argument Parser (argparse)
    parser = argparse.ArgumentParser(description="Generate and train a Neural Network")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--verbose", action="store_true",
                       help="Increase Verbosity.")
    group.add_argument("-q", "--quiet", action="store_true",
                       help="Decrease Verbosity.")
    parser.add_argument("--no-confirm", action="store_false",
                        help="Use to skip confirmation before training.")
    parser.add_argument("-cf", "--config-file", default=None, type=str,
                        help="Path to config file")
    parser.add_argument("-rcf", "--reset-config-file", action="store_true",
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
    args = parser.parse_args()
    if args.reset_config_file and args.config_file is None:
        parser.error("--reset_config_file requires --config_file.")

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
    # data_path = parameters["data_path"]
    data_path = ""
    data_sufix = parameters["data_sufix"]
    # output_model_path = parameters["output_model_path"]
    # path_output_plots = "/Volumes/Elements/Outputs/plot"
    # output_plots_path = "plot"

    # Load Data
    print("[INFO] Loading data...")
    input_spectra, input_magnitudes, label_sfh, label_z, spectra_lambda, agevec = \
        loadfiles(input_path=data_path + "Input_" + data_sufix + ".fits",
                  labels_path=data_path + "Label_" + data_sufix + ".fits")

    print(f"""
    Variable sizes:
        Input_spectra: {input_spectra.shape} - {convert_bytes(input_spectra.nbytes)}
        Input_magnitudes: {input_magnitudes.shape} - {convert_bytes(input_magnitudes.nbytes)}
        Label_sfh: {label_sfh.shape} - {convert_bytes(label_sfh.nbytes)}
        Label_z: {label_z.shape} - {convert_bytes(label_z.nbytes)}
        """)

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
    train_data = np.concatenate([trainSpect, trainMag], axis=1)
    train_labels = np.concatenate([trainLabSfh, trainLabZ], axis=1)
    test_data = np.concatenate([testSpect, testMag], axis=1)
    test_labels = np.concatenate([testLabSfh, testLabZ], axis=1)

    # Verify the size of the
    print_train_test_sizes(split_train_test, main_title="Sizes of diferent sets (Train, Test)")
    # "split" tuple can be removed now
    del split_train_test

    ############################################################################
    # Build model
    # custom_kwargs_model = {"spect_neurons_first_layer": 256}
    # model = Cerebro.build_model(epochs=epochs, loss_function_used=loss_function_used, init_lr=init_lr,
    #                             **custom_kwargs_model, **main_kwargs)
    # model.summary()
    # Cerebro.graph(model, "tstimage4.png")

    ############################################################################
    # Cross Validation
    # Define Regressor Model

    # Set parameters that will generate the grid
    # Here is where all the parameters will go (regarding 'spect_', ...)
    # param_grid = dict(
    #     epochs=[3, 4],
    #     batch_size=[200, 300],
    #     magn_neurons_first_layer=[32, 64],
    #     magn_number_layers=[2, 3],
    #     spect_number_layers=[4],
    #     spect_neurons_first_layer=[128, 256],
    #     spect_filter_size=[[30, 20, 10, 5], [30, 30, 10, 5], [50, 25, 10, 5], [30, 15, 10, 5], 15],
    # )
    param_grid = cv_parameters
    print("param_grid", param_grid)
    number_of_combinations = 1
    for key, value in param_grid.items():
        print(f"{key}: {value} - {len(value)}")
        number_of_combinations *= len(value)
    print(f"[INFO] Number of possible combinations: {number_of_combinations}")

    estimator = KerasRegressor(model=Cerebro.build_model, **parameters, **param_grid)
    print("estimator", estimator.get_params().keys())

    grid = RandomizedSearchCV(estimator=estimator, param_distributions=param_grid,
                              n_iter=2,
                              cv=cv, random_state=traintestrandomstate, verbose=2)
    # grid = GridSearchCV(estimator=estimator, param_grid=param_grid, cv=cv, verbose=5)
    # print("grid", grid.get_params().keys())

    # Train model
    if not parameters["no_confirm"]:
        continue_with_training = input("Continue training? [Y]/N ")
        if continue_with_training.lower() in ["n", "no", "stop"]:
            raise KeyboardInterrupt('User stopped the execution.')
    print("[INFO] Start training...")
    grid_result = grid.fit(train_data, train_labels)













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

    accuracy = grid.score(test_data, test_labels)
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
    main(config_file_path="config.txt")
