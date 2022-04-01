# -*- coding: utf-8 -*-
# Just A Network Executor [Phoenix]


# help("modules")
# help("modules tensorflow")
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import random
from sklearn.model_selection import train_test_split
from sklearn.model_selection import RandomizedSearchCV
from scikeras.wrappers import KerasRegressor
import joblib
import keras
import sys
from pprint import pprint
from keras.callbacks import EarlyStopping, ModelCheckpoint

from MUTANTS.ERIK import *
from MUTANTS.RAVEN import convert_bytes, print_train_test_sizes
from MUTANTS.XAVIER import Cerebro

print("[INFO] Finished Importing")


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
    # ToDo: Fix this in the final version (i.e. no more debug)
    # data_path = parameters["data_path"]
    # data_path = ""
    data_path = "/Volumes/Elements/Outputs/"
    # parameters["output_model_path"] = ""
    # parameters["output_model_path"] = "/Volumes/Elements/Outputs/"
    # data_sufix = parameters["data_sufix"]
    data_sufix = "combined"
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
    Cerebro.graph(model, "tstimage.png")

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
            ax[i].plot(np.arange(0, epochs_draw), history.history[l], label=l)
            ax[i].plot(np.arange(0, epochs_draw), history.history["val_" + l], label="val_" + l)
            ax[i].legend()
        # # save the losses figure
        # plt.tight_layout()
        # plt.savefig("{}_losses.png".format(path_output_plots))
        # print("[INFO] Loss image stored in {}_losses.png".format(path_output_plots))
        # plt.close()

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



    # load the saved model
    saved_model = keras.models.load_model(best_model_path, custom_objects={"smape_loss": Cerebro.smape_loss})
    # evaluate the model
    print("[Evaluate] Evaluating best model")
    print("skipping")
    # _, train_acc = saved_model.evaluate(train_data, train_labels, verbose=1)
    # _, test_acc = saved_model.evaluate(test_data, test_labels, verbose=1)
    # print('Train: %.3f, Test: %.3f' % (train_acc, test_acc))










    plt.show()

    sys.exit(2)
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
    print("param_grid", param_grid, "\n")
    number_of_combinations = 1
    for key, value in param_grid.items():
        print(f"{key}: {value} - {len(value)}")
        number_of_combinations *= len(value)
    print(f"[INFO] Number of possible combinations: {number_of_combinations}")

    estimator = KerasRegressor(model=Cerebro.build_model, **parameters, **param_grid)
    # ToDo: Remove these lines
    print("\n[INFO] THIS IS EVERYTHING THAT IS GOING TO GO INSIDE, THE PREVIOUS LINES ARE THE CORRECT ONES")
    print("estimator", estimator.get_params().keys())

    grid = RandomizedSearchCV(estimator=estimator, param_distributions=param_grid,
                              n_iter=parameters["n_iter"],
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
    try:
        print(cv_results["mean_test_score"][1] - cv_results["mean_test_score"][2])
    except Exception as ex:
        template = "An exception of type {0} occurred. Arguments:\n{1!r}"
        message = template.format(type(ex).__name__, ex.args)
        print(message)
        print("[ERROR] If the error is of type KeyError, the issue was that there are too few cases cv cases!")

    # Saving Model
    print(f"[INFO] Saving fitted object to folder '{parameters['output_model_path']}'")
    print("[INFO] Saving grid_result to", parameters["output_name_grid_result"], "...")
    joblib.dump(grid_result, os.path.join(parameters["output_model_path"],
                                          parameters["output_name_grid_result"]))

    print("[INFO] Saving best estimator to", parameters["output_name_best_estimator"], "...")
    joblib.dump(grid_result.best_estimator_, os.path.join(parameters["output_model_path"],
                                                          parameters["output_name_best_estimator"]), compress=1)

    print("[INFO] Saving best model to", parameters["output_name_model"], "...")
    grid.best_estimator_.model_.save(os.path.join(parameters["output_model_path"],
                                                  parameters["output_name_model"]))

    print("[INFO] Saving parameters to", parameters["output_name_parameters"], "...")
    output_params = {"parameters": parameters, "param_grid": param_grid, "best_params": grid_result.best_params_}
    json.dump(output_params, open(os.path.join(parameters["output_model_path"],
                                               parameters["output_name_parameters"]), "w+"))

    # Load Model for testing
    # IMPORTANT: when loading, the custom loss function needs to be given!
    model_loaded = keras.models.load_model(os.path.join(parameters["output_model_path"],
                                                        parameters["output_name_model"]),
                                           custom_objects={"smape_loss": Cerebro.smape_loss})
    parameters_loaded = json.load(open(os.path.join(parameters["output_model_path"],
                                                    parameters["output_name_parameters"]), "r"))
    print("Evaluate on test data")
    results = model_loaded.evaluate(test_data, test_labels, batch_size=parameters_loaded["best_params"]["batch_size"])
    print("test loss, test acc:", results)

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
    # print(getparametersfromid("/Volumes/Elements/Outputs/MetaD_combined.json", 34509, verbose=0))
    # combine_datasets(["20220222T112103_drDhe7", "20220222T134553_qbkxpF", "20220222T213001_nHS1Bf"],
    #                  "/Volumes/Elements/Outputs/",
    #                  )
    main(config_file_path="config.txt")
