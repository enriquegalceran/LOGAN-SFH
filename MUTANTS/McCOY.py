# -*- coding: utf-8 -*-
# Measurements, COrrelation and sanitY of data

import argparse
import os
import sys
import csv
from typing import Callable, Any

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import keras
import tensorflow
from sklearn.model_selection import train_test_split

import ERIK
import XAVIER
import RAVEN


def cleanlogfile(filename, filename_out=None):
    with open(filename, 'r') as f:
        lines = f.readlines()

    lines = [line.strip() for line in lines if "] - ETA:" not in line]

    for line in lines:
        print(line)

    lines2 = [i for i in lines if "=] - " in i]
    lines2 = [i.split("] - ")[1] for i in lines2]
    lines2 = [int(i.split("s ", 1)[0]) for i in lines2]

    print("times per epoch:")
    print(lines2)
    print("sum:", sum(lines2))
    print("min:", sum(lines2) / 60)
    print("hours:", sum(lines2) / 3600)

    lines.append(f"s_total: {sum(lines2)}")
    lines.append(f"m_total: {sum(lines2) / 60}")
    lines.append(f"h_total: {sum(lines2) / 3600}")

    lines = [line + "\n" for line in lines]
    with open(filename_out, "w+") as f:
        f.writelines(lines)


def model_colormap(model_paths, data_path, name_models=None, loss=None,
                   method_standardize_label_sfh=0,
                   method_standardize_label_z=0,
                   method_standardize_spectra=0,
                   method_standardize_magnitudes=0,
                   which_data_is_going_to_be_used=2,
                   first_id_plot=0, n_plot=None,
                   mageburst=1e8,
                   train_size=0.8, test_size=0.2, traintestrandomstate=42, traintestshuffle=True):

    # Load Data
    print("[INFO] Loading data...")
    label_path = data_path.replace("Input_", "Label_")
    metadata_path = data_path.replace("Input_", "MetaD_").replace(".fits", ".rda")

    input_spectra, input_magnitudes, label_sfh, label_z, spectra_lambda, agevec, ageweight = \
        ERIK.loadfiles(input_path=data_path, labels_path=label_path,
                       method_standardize_label_sfh=method_standardize_label_sfh,
                       method_standardize_label_z=method_standardize_label_z,
                       method_standardize_spectra=method_standardize_spectra,
                       method_standardize_magnitudes=method_standardize_magnitudes)

    _, _, label_sfh_real, label_z_real, _, _, _ = \
        ERIK.loadfiles(input_path=data_path, labels_path=label_path,
                       method_standardize_label_sfh=0,
                       method_standardize_label_z=0,
                       method_standardize_spectra=0,
                       method_standardize_magnitudes=0)

    # Verified: If indices are added, they will STILL respect the output as if there were no indices
    split_train_test = train_test_split(input_spectra, input_magnitudes, label_sfh, label_z, label_sfh_real,
                                        label_z_real,
                                        range(input_spectra.shape[0]),  # Indices
                                        test_size=test_size,
                                        train_size=train_size,
                                        random_state=traintestrandomstate,
                                        shuffle=traintestshuffle)
    (trainSpect, testSpect,
     trainMag, testMag,
     trainLabSfh, testLabSfh,
     trainLabZ, testLabZ,
     trainLabSfh_real, testLabSfh_real,
     trainLabZ_real, testLabZ_real,
     trainIndices, testIndices) = split_train_test

    if name_models is None:
        name_models = model_paths

    saved_models = []
    for pqwe, model_p in enumerate(model_paths):
        if pqwe == 1:
            break
        saved_models.append(keras.models.load_model(model_p,
                                                    custom_objects={"smape_loss": XAVIER.Cerebro.smape_loss}))
        # saved_models[0].summary()

    #################
    # Decide which data is going to be used
    if which_data_is_going_to_be_used == 0:
        # All the data
        in_data = np.concatenate([input_spectra, input_magnitudes], axis=1)
        idx = list(range(input_spectra.shape[0]))
    elif which_data_is_going_to_be_used == 1:
        # Only Training data
        in_data = np.concatenate([trainSpect, trainMag], axis=1)
        label_sfh = trainLabSfh
        label_z = trainLabZ
        label_sfh_real = trainLabSfh_real
        label_z_real = trainLabZ_real
        idx = trainIndices
    elif which_data_is_going_to_be_used == 2:
        # Only Test data
        in_data = np.concatenate([testSpect, testMag], axis=1)
        label_sfh = testLabSfh
        label_z = testLabZ
        label_sfh_real = testLabSfh_real
        label_z_real = testLabZ_real
        idx = testIndices
    else:
        raise ValueError(f"which_data_is_going_to_be_used should be [0, 1, 2] and is {which_data_is_going_to_be_used}")

    if loss is None:
        loss = tensorflow.keras.losses.MeanSquaredError()

    if n_plot is None:
        n_plot = in_data.shape[0] - first_id_plot

    outputs = []
    losses = []
    for model in saved_models:
        outputs.append(model.predict(in_data[first_id_plot:(first_id_plot + n_plot)]))
        losses.append(
            [model.evaluate(np.array([in_data[idx, :]]), [np.array([label_sfh[idx, :]]), np.array([label_z[idx, :]])],
                            batch_size=1, verbose=1)
             for idx in range(first_id_plot, first_id_plot + n_plot)]
        )

    combined_losses = [[k[0] for k in losses[i_model]] for i_model in range(len(saved_models))]
    log_combined_losses = [np.log10(_) for _ in combined_losses]

    metallicity_loss = [[k[1] for k in losses[i_model]] for i_model in range(len(saved_models))]
    log_metallicity_losses = [np.log10(_) for _ in metallicity_loss]

    sfh_losses = [[k[2] for k in losses[i_model]] for i_model in range(len(saved_models))]
    log_sfh_losses = [np.log10(_) for _ in sfh_losses]

    # mass_weighted = label_sfh * ageweight
    # total_mass = np.sum(mass_weighted, axis=1)
    mass_weigted = label_sfh_real * ageweight
    total_mass_no_norm = np.sum(mass_weigted, axis=1)
    mass_only_burst_idx = agevec < mageburst
    mass_burst = np.sum(mass_weigted[:, mass_only_burst_idx], axis=1)

    def plot_colorbars(x, y, c, s=3, cmap='plasma_r', log="xy", xlim=None, ylim=None,
                       xlabel=None, ylabel=None, collabel=None, title=None,
                       savefig_name=None, figsize=(12, 9)):
        fig, ax = plt.subplots(figsize=figsize)
        plt.title(title)
        sc = ax.scatter(x, y, c=c, s=s, cmap=cmap)
        if 'x' in log:
            ax.set_yscale('log')
        if "y" in log:
            ax.set_xscale('log')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        cbar = plt.colorbar(sc, ax=ax)
        cbar.set_label(collabel)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        plt.tight_layout()
        if savefig_name is not None:
            plt.savefig(savefig_name)
        plt.show()

    plot_colorbars(x=total_mass_no_norm[first_id_plot:(first_id_plot + n_plot)],
                   y=mass_burst[first_id_plot:(first_id_plot + n_plot)],
                   c=log_combined_losses[0],
                   xlabel=r"$M_{total} (M_{\odot})$",
                   ylabel=r"$M_{burst} (M_{\odot}) (\leq0.1Gyr)$",
                   collabel=r'$Log_{10}(\ell)$',
                   title=r"Combined Loss Function ($\ell$) based on $M_{total}$ and $M_{burst}$",
                   savefig_name='/Users/enrique/Documents/GitHub/LOGAN-SFH/my_test_fig_combined.png',
                   )
    plot_colorbars(x=total_mass_no_norm[first_id_plot:(first_id_plot + n_plot)],
                   y=mass_burst[first_id_plot:(first_id_plot + n_plot)],
                   c=log_sfh_losses[0],
                   xlabel=r"$M_{total} (M_{\odot})$",
                   ylabel=r"$M_{burst} (M_{\odot}) (\leq0.1Gyr)$",
                   collabel=r'$Log_{10}(\ell)$',
                   title=r"SFH Loss Function ($\ell$) based on $M_{total}$ and $M_{burst}$",
                   savefig_name='/Users/enrique/Documents/GitHub/LOGAN-SFH/my_test_fig_sfh.png',
                   )
    plot_colorbars(x=total_mass_no_norm[first_id_plot:(first_id_plot + n_plot)],
                   y=mass_burst[first_id_plot:(first_id_plot + n_plot)],
                   c=log_metallicity_losses[0],
                   xlabel=r"$M_{total} (M_{\odot})$",
                   ylabel=r"$M_{burst} (M_{\odot}) (\leq0.1Gyr)$",
                   collabel=r'$Log_{10}(\ell)$',
                   title=r"Metallicity Loss Function ($\ell$) based on $M_{total}$ and $M_{burst}$",
                   savefig_name='/Users/enrique/Documents/GitHub/LOGAN-SFH/my_test_fig_metal.png',
                   )

   # xlim=(50327.76943032607, 368356365924.79095),
   # ylim=(60.62030449954481, 819673467.9265951))

    print("------")















def verify_model(model_paths, data_path, name_models=None, loss=None,
                 method_standardize_label_sfh=0,
                 method_standardize_label_z=0,
                 method_standardize_spectra=0,
                 method_standardize_magnitudes=0,
                 which_data_is_going_to_be_used=2,
                 first_id_plot=0, n_plot=3,
                 train_size=0.8, test_size=0.2, traintestrandomstate=42, traintestshuffle=True):
    # Load Data
    print("[INFO] Loading data...")
    label_path = data_path.replace("Input_", "Label_")
    metadata_path = data_path.replace("Input_", "MetaD_").replace(".fits", ".rda")

    input_spectra, input_magnitudes, label_sfh, label_z, spectra_lambda, agevec, _ = \
        ERIK.loadfiles(input_path=data_path, labels_path=label_path,
                       method_standardize_label_sfh=method_standardize_label_sfh,
                       method_standardize_label_z=method_standardize_label_z,
                       method_standardize_spectra=method_standardize_spectra,
                       method_standardize_magnitudes=method_standardize_magnitudes)

    _, _, label_sfh_real, label_z_real, _, _, _ = \
        ERIK.loadfiles(input_path=data_path, labels_path=label_path,
                       method_standardize_label_sfh=0,
                       method_standardize_label_z=0,
                       method_standardize_spectra=0,
                       method_standardize_magnitudes=0)

    # Verified: If indices are added, they will STILL respect the output as if there were no indices
    split_train_test = train_test_split(input_spectra, input_magnitudes, label_sfh, label_z, label_sfh_real, label_z_real,
                                        range(input_spectra.shape[0]),      # Indices
                                        test_size=test_size,
                                        train_size=train_size,
                                        random_state=traintestrandomstate,
                                        shuffle=traintestshuffle)
    (trainSpect, testSpect,
     trainMag, testMag,
     trainLabSfh, testLabSfh,
     trainLabZ, testLabZ,
     trainLabSfh_real, testLabSfh_real,
     trainLabZ_real, testLabZ_real,
     trainIndices, testIndices) = split_train_test

    if name_models is None:
        name_models = model_paths

    saved_models = []
    for model_p in model_paths:
        saved_models.append(keras.models.load_model(model_p,
                                                    custom_objects={"smape_loss": XAVIER.Cerebro.smape_loss}))
        # saved_models[0].summary()

    #################
    # Decide which data is going to be used
    if which_data_is_going_to_be_used == 0:
        # All the data
        in_data = np.concatenate([input_spectra, input_magnitudes], axis=1)
        idx = list(range(input_spectra.shape[0]))
    elif which_data_is_going_to_be_used == 1:
        # Only Training data
        in_data = np.concatenate([trainSpect, trainMag], axis=1)
        label_sfh = trainLabSfh
        label_z = trainLabZ
        label_sfh_real = trainLabSfh_real
        label_z_real = trainLabZ_real
        idx = trainIndices
    elif which_data_is_going_to_be_used == 2:
        # Only Test data
        in_data = np.concatenate([testSpect, testMag], axis=1)
        label_sfh = testLabSfh
        label_z = testLabZ
        label_sfh_real = testLabSfh_real
        label_z_real = testLabZ_real
        idx = testIndices
    else:
        raise ValueError(f"which_data_is_going_to_be_used should be [0, 1, 2] and is {which_data_is_going_to_be_used}")

    if loss is None:
        loss = tensorflow.keras.losses.MeanSquaredError()

    outputs = []
    for model in saved_models:
        outputs.append(model.predict(in_data[first_id_plot:(first_id_plot + n_plot)]))

    legend_name = ["Label"] + name_models
    for i in range(first_id_plot, first_id_plot + n_plot):
        fig, ax = plt.subplots(2, 2, figsize=(25, 15))
        plt.suptitle(f"i:{i + 1} - ID:{idx[i]}")
        ax[0, 0].set_title("SFH")
        ax[0, 0].plot(agevec, label_sfh[i, :], 'k')
        for q, k in enumerate(outputs):
            ax[0, 0].plot(agevec, k[0][(i - first_id_plot), :])

        ax[0, 0].set_xscale('log')
        ax[0, 0].legend(legend_name)

        ax[1, 0].set_title("SFH - residuals")
        ax[1, 0].plot(agevec, np.zeros(agevec.shape), 'k')
        current_legend = ["Label"]
        for q, k in enumerate(outputs):
            ax[1, 0].scatter(agevec, np.subtract(k[0][(i - first_id_plot), :], label_sfh[i, :]))
            current_legend.append(name_models[q] + "_" +
                                  f"{loss(label_sfh[i, :], k[0][(i - first_id_plot), :]).numpy():3f}")
        ax[1, 0].set_xscale('log')
        ax[1, 0].legend(current_legend)

        ax[0, 1].set_title("Metallicity")
        ax[0, 1].plot(agevec, label_z[i, :], 'k')
        for k in outputs:
            ax[0, 1].plot(agevec, k[1][(i - first_id_plot), :])
        ax[0, 1].set_xscale('log')
        ax[0, 1].legend(legend_name)

        ax[1, 1].set_title("Metallicity - residuals")
        ax[1, 1].plot(agevec, np.zeros(agevec.shape), 'k')
        current_legend = ["Label"]
        for q, k in enumerate(outputs):
            ax[1, 1].scatter(agevec, np.subtract(k[1][(i - first_id_plot), :], label_z[i, :]))
            current_legend.append(name_models[q] + "_" +
                                  f"{loss(label_z[i, :], k[1][(i - first_id_plot), :]).numpy():4f}")
        ax[1, 1].set_xscale('log')
        ax[1, 1].legend(current_legend)
        plt.tight_layout()
        plt.show()

        tmp_dic = {"agevec": agevec,
                   "sfh_true": label_sfh[i, :],
                   "z_true": label_z[i, :],
                   "sfh_no_stand": label_sfh_real[i, :],
                   "z_no_stand": label_z_real[i, :]}
        for q, k in enumerate(outputs):
            tmp_dic[f"sfh{q}"] = k[0][(i - first_id_plot), :]
            tmp_dic[f"z{q}"] = k[1][(i - first_id_plot), :]
        tmp_sp = {"spectr_in": input_spectra[i, :],
                  "waveout": spectra_lambda}
        tmp_mag = {"magnitudes_in": input_magnitudes[i, :],
                   "ID": [idx[i]]*5}
        tmp_id_names = {"ID": [idx[i]]}
        for q, k in enumerate(legend_name):
            if q == 0:
                continue
            else:
                tmp_id_names[f"name{q-1}"] = [k]
        tmp_df = pd.DataFrame(tmp_dic)
        tmp_df.to_csv("/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file.pd", index=False)
        tmp_df2 = pd.DataFrame(tmp_sp)
        tmp_df3 = pd.DataFrame(tmp_mag)
        tmp_df4 = pd.DataFrame(tmp_id_names)
        tmp_df2.to_csv("/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file2.pd", index=False)
        tmp_df3.to_csv("/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file3.pd", index=False)
        tmp_df4.to_csv("/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file4.pd", index=False)





        print(i)

    return
    difference_sfh = np.subtract(output[0], label_sfh[first_id_plot:first_id_plot + n_plot, :])
    difference_mag = np.subtract(output[1], label_z[first_id_plot:first_id_plot + n_plot, :])

    fig1, ax1 = plt.subplots()
    ax1.set_title("Boxplot SFH")
    ax1.boxplot(difference_sfh)
    fig2, ax2 = plt.subplots()
    ax2.set_title("Boxplot Z")
    ax2.boxplot(difference_mag)
    fig1, ax1 = plt.subplots()
    ax1.set_title("Boxplot SFH (relativo)")
    ax1.boxplot(np.divide(difference_sfh, label_sfh[first_id_plot:first_id_plot + n_plot, :]))
    fig2, ax2 = plt.subplots()
    ax2.set_title("Boxplot Z (relativo)")
    ax2.boxplot(np.divide(difference_mag, label_z[first_id_plot:first_id_plot + n_plot, :]))

    # for ids in range(10):
    #
    #     metadata = ERIK.getparametersfromid("/Users/enrique/Documents/GitHub/LOGAN-SFH/MetadataOutput.json",
    #                                         ids, verbose=0, returnfunction=False)
    #     fig3, ax3 = plt.subplots()
    #     ax3.set_title(f"SFH - {ids}")
    #     ax3.semilogx(agevec, label_sfh[5000 + ids], "k-")
    #     ax3.semilogx(agevec, output[0][ids], "g.-")
    #     plt.savefig(f"/Users/enrique/Documents/GitHub/LOGAN-SFH/Test/SFH{ids}.png")
    #     fig4, ax4 = plt.subplots()
    #     ax4.set_title(f"Z - {ids}")
    #     ax4.semilogx(agevec, label_z[5000 + ids], "k-")
    #     ax4.semilogx(agevec, output[1][ids], "g.-")
    #     plt.savefig(f"/Users/enrique/Documents/GitHub/LOGAN-SFH/Test/Z{ids}.png")

    plt.show()
    print("test1234")


def main(update_values=None, **mainkwargs):
    # Argument Parser (argparse)
    parser = argparse.ArgumentParser(description="Cleaning programm")

    parser.add_argument("-ln", "--logname", default=None, type=str,
                        help="If a filename is given, it will clean this log (i.e. remove extra characters")
    parser.add_argument("-clo", "--cleanlogoutput", default=None, type=str,
                        help="If a filename is given, will write the cleaned log in this location. "
                             "Requires a value in --logname to work.")
    parser.add_argument("-vm", "--verify_model", default=None, type=list, nargs="*",
                        help="Test the prediction of the Neural Network. the absolute/relative path to the model is"
                             " required. Requieres a dataset to verify (--data_path).")
    parser.add_argument("-dp", "--data_path", default=None, type=str,
                        help="Indicates a dataset to verify the model. If the data is not a path (absolute/relative) "
                             "(if it does not end in '.fits'), it will consider the data to be located in the current "
                             "folder and that the given value is the sufix. REQUIRES THE INPUT DATA. The Labels are "
                             "expected to be located in the same folder.")
    parser.add_argument("-v", "--verbose", default=1, type=int,
                        help="Indicates the value of verbosity. Minimum=0, Default=1.")
    parser.add_argument("-jc", "--json_clean", action="append", nargs="+", default=None, type=str,
                        help="Path to json file to be cleaned. Multiple files can be given (either using -jc again, or "
                             "be simply giving multiple paths).If verbose > 0, json file will be printed.")
    args = parser.parse_args()

    # Set parameters passed through main
    if update_values is not None:
        for key, value in update_values.items():
            args.__setattr__(key, value)

    # Verify parameters
    if args.cleanlogoutput is not None and args.logname is None:
        parser.error("--cleanlogoutput requires --logname")
    if args.verify_model is not None and args.data_path is None:
        parser.error("--verify_model requieres --data_path")
    if args.verify_model is None and args.data_path is not None:
        parser.error("--data_path requieres --verify_model")

    if args.logname is not None:
        print(f"[INFO] Cleaning log '{args.logname}'")
        cleanlogfile(args.logname, args.cleanlogoutput)
        print("[LOG] FINISHED")

    if args.verify_model is not None:
        # If only one element, convert to list
        if type(args.verify_model) is str:
            args.verify_model = [args.verify_model]
        print(f"[INFO] Verifying Model(s) {args.verify_model} with Data: {args.data_path}")
        data_path = args.data_path
        if data_path[-5:] != ".fits":
            data_path = "Input_" + data_path + ".fits"
        if data_path[0] != "/":
            data_path = os.path.abspath(data_path)

        for i, model in enumerate(args.verify_model):
            if model[0] != "/":
                args.verify_model[i] = os.path.abspath(model)

        # Verify Model
        model_names = None
        if "model_names" in mainkwargs:
            model_names = mainkwargs["model_names"]
        arguments_standardize = {}
        for key in mainkwargs.keys():
            if "method_standardize_" in key:
                arguments_standardize[key] = mainkwargs[key]
        model_colormap(args.verify_model, data_path, model_names, **arguments_standardize)
        # verify_model(args.verify_model, data_path, model_names, **arguments_standardize)

    if args.json_clean is not None:
        json_files = RAVEN.flatten_list(args.json_clean)
        for file in json_files:
            ERIK.prettyfy_json_file(file, args.verbose)


if __name__ == "__main__":
    # models_path = "/Users/enrique/Documents/GitHub/LOGAN-SFH/TrainedModels/MSE_reduced_wave_small_dataset/"
    models_path = "/Users/enrique/scpdata/2"
    data_path = os.path.join("/Users/enrique/Documents/GitHub/LOGAN-SFH/TrainingData/",
                             "Input_combined_wave.fits")

    files = os.listdir(models_path)
    keywords_to_remove = ["epoch001", "imagen", "config"]
    for keyword in keywords_to_remove:
        files = [_ for _ in files if keyword not in _]
    files.sort()
    length_prefix_epoch = files[0].find("epoch")
    length_prefix_loss = files[0].find("loss")
    filenames = [f"e{_[(length_prefix_epoch + 5):(length_prefix_epoch + 8)]}-"
                 f"l{_[(length_prefix_loss + 4):-3]}" for _ in files]
    models = [os.path.join(models_path, f) for f in files]

    main({"verify_model": models, "data_path": data_path},
         model_names=filenames,
         method_standardize_label_sfh=5,
         method_standardize_label_z=0
         )

    # old_models = ["/Users/enrique/model_mean_squared_13_05.h5",
    #               "/Users/enrique/model_mean_squared_log_13_05.h5",
    #               "/Users/enrique/Documents/GitHub/LOGAN-SFH/model_no_normalizing_12_05.h5"]

    # main({"verify_model": os.path.join("/Users/enrique/Documents/GitHub/LOGAN-SFH/",
    #                                    "model_no_normalizing_12_05.h5"),
    #       "data_path": os.path.join("/Users/enrique/Documents/GitHub/LOGAN-SFH/TrainingData/",
    #                                 "Input_combined.fits")})


# ToDo: comparar dos SFH (ancha y estrecha), sacar gráfica de SFH y metalicidades
# Mandar el código de la red
# meter la salida de la red en prospect y mirar el espectro+magnitudes que obtenemos.
# Probar SIN ruido?
