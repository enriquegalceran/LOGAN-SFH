# -*- coding: utf-8 -*-
# Measurements, COrrelation and sanitY of data

import argparse
import os
# import sys
# import csv
# from typing import Callable, Any
from pylick.pylick.analysis import Galaxy
from pylick.pylick.indices import IndexLibrary
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import keras
# import tensorflow
from sklearn.model_selection import train_test_split
from math import nan

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


def save_data_to_file(model_paths, data_path,
                      method_standardize=None,
                      which_data="val",
                      first_id_plot=0, n_plot=None,
                      temp_folder="/Users/enrique/Documents/GitHub/LOGAN-SFH/tempFolder",
                      **kwargs):
    # Plan: export data (csv?) for R to be read
    # Load Data
    loaded_data = evaluate_model(model_paths, data_path, return_loss=False, method_standardize=method_standardize,
                                 which_data=which_data, first_id_plot=first_id_plot, n_plot=n_plot, **kwargs)
    (input_spectra, input_magnitudes, label_sfh, label_z, label_sfh_no_normalization, label_z_no_normalization,
     spectra_lambda, agevec, ageweight, saved_models, idx, outputs, losses) = loaded_data
    # remove tuple data
    del loaded_data

    print(f"[INFO] Exported files to folder: {os.path.abspath(temp_folder)}")
    # Export all the relevant information in numpy format in the temp_folder directory
    np.save(os.path.join(temp_folder, "input_spectra"), input_spectra)
    np.save(os.path.join(temp_folder, "input_magnitudes"), input_magnitudes)
    np.save(os.path.join(temp_folder, "agevec"), agevec)
    np.save(os.path.join(temp_folder, "ageweight"), ageweight)
    np.save(os.path.join(temp_folder, "wave"), spectra_lambda)
    np.save(os.path.join(temp_folder, "label_sfh"), label_sfh)
    np.save(os.path.join(temp_folder, "label_z"), label_z)
    np.save(os.path.join(temp_folder, "label_sfh_no_normalization"), label_sfh_no_normalization)
    np.save(os.path.join(temp_folder, "label_z_no_normalization"), label_z_no_normalization)
    for model_id in range(len(saved_models)):
        np.save(os.path.join(temp_folder, f"predict_sfh_{model_id}"), outputs[model_id][0])
        np.save(os.path.join(temp_folder, f"predict_z_{model_id}"), outputs[model_id][1])
    # ToDo: Add losses here? They need a bit more of preprocessing before being saved, so I left this for afterwards.


def calculate_lick_indices_mass(idx_models=None, calculate_for_input=True, force_error_sn=100,
                                temp_folder="/Users/enrique/Documents/GitHub/LOGAN-SFH/tempFolder",
                                return_indices=True, calculate_differences=True, verbosity_progress=100,
                                **kwargs):
    # ToDo: Acknoledgements need to be added: https://gitlab.com/mmoresco/pylick#index_list
    if idx_models is None:
        idx_models = [0, 1]

    def single_relative_difference(t, p):
        return (t - p) / t

    def list_relative_difference(true, predicted):
        return [single_relative_difference(true[k], predicted[k]) for k in range(len(true))]

    # Name of files to be loaded
    name_input_spectra = os.path.join(temp_folder, "input_spectra.npy")
    name_wave = os.path.join(temp_folder, "wave.npy")
    name_files_predicted = [os.path.join(temp_folder, f"predict_spectra_{id_mod}.npy") for id_mod in idx_models]
    index_values_predicted = None
    index_values_spectra = None
    predicted_spectra_number = None
    out = {}

    # indexes that are going to be used:
    index_list = np.arange(53, 66)

    # Load wave
    wave = np.load(name_wave)

    # Calculate index for input_spectra
    if calculate_for_input:
        input_spectra = np.load(name_input_spectra)
        input_spectra_number = input_spectra.shape[0]
        index_values_spectra = np.zeros((input_spectra_number, len(index_list)), dtype=float)

        print("[INFO] Calculating lick indices for input spectra...")
        start_time = datetime.now()
        for i in range(input_spectra_number):
            if i % verbosity_progress == 0:
                print(f"Current case for input spectra : {i}/{input_spectra_number}")
            tmp = Galaxy("Input_" + str(i), index_list, spec_wave=wave, spec_flux=input_spectra[i, :],
                         spec_err=input_spectra[i, :] / force_error_sn, spec_mask=None, meas_method='int')
            index_values_spectra[i, :] = tmp.vals
        end_time = datetime.now()
        print('Duration: {}'.format(end_time - start_time))
        print("[INFO] Finished calculating lick indices for input spectra.")
        print(f"[INFO] storing in {os.path.join(temp_folder, 'lick_idx_input_spectra.npy')} ...")
        np.save(os.path.join(temp_folder, "lick_idx_input_spectra.npy"), index_values_spectra)
    else:
        index_values_spectra = np.load(os.path.join(temp_folder, "lick_idx_input_spectra.npy"))

    for imodel in idx_models:
        predicted_spectra = np.load(name_files_predicted[imodel])
        predicted_spectra_number = predicted_spectra.shape[0]
        index_values_predicted = np.zeros((predicted_spectra_number, len(index_list)), dtype=float)

        print(f"[INFO] Calculating lick indices for predicted spectra for model #{imodel}...")
        start_time = datetime.now()
        for i in range(predicted_spectra_number):
            if i % verbosity_progress == 0:
                print(f"Current case for predicted spectra model #{imodel} : {i}/{predicted_spectra_number}")
            # FIXME: There are some cases where there are NaNs (maybe because sfh/z is predicted as <0?)
            if not np.isnan(predicted_spectra[i, 0]):
                tmp = Galaxy("Predicted_" + str(i), index_list, spec_wave=wave, spec_flux=predicted_spectra[i, :],
                             spec_err=predicted_spectra[i, :] / force_error_sn, spec_mask=None, meas_method='int')
                index_values_predicted[i, :] = tmp.vals
            else:
                index_values_predicted[i, :] = [nan] * len(index_list)
        end_time = datetime.now()
        print('Duration: {}'.format(end_time - start_time))
        print(f"[INFO] Finished calculating lick indices for predicted spectra for model #{imodel}.")
        print(f"[INFO] storing in {os.path.join(temp_folder, f'lick_idx_predicted_{imodel}.npy')} ...")
        np.save(os.path.join(temp_folder, f"lick_idx_predicted_{imodel}.npy"), index_values_predicted)

        if calculate_differences:
            differences = np.zeros((predicted_spectra_number, ))
            print(f"[INFO] Calculating differences for predicted spectra vs. input spectra for model #{imodel}...")
            start_time = datetime.now()
            for i in range(predicted_spectra_number):
                if i % verbosity_progress == 0:
                    print("Current case for differences for predicted spectra vs. input spectra model "
                          f"#{imodel} : {i}/{predicted_spectra_number}")
                differences[i] = np.square(list_relative_difference(index_values_spectra[i, ],
                                                                    index_values_predicted[i, ])
                                           ).mean()
            end_time = datetime.now()
            print('Duration: {}'.format(end_time - start_time))
            print("[INFO] Finished calculating differences for predicted spectra vs. input "
                  f"spectra for model #{imodel}.")
            print(f"[INFO] storing in {os.path.join(temp_folder, f'differences_{imodel}.npy')} ...")
            np.save(os.path.join(temp_folder, f"differences_{imodel}.npy"), differences)
            out[f"Relative-MSE_{imodel}"] = differences

    if return_indices:
        out["index_values_predicted"] = index_values_predicted
        out["index_values_spectra"] = index_values_spectra

    return out


def evaluate_model(model_paths, data_path, return_loss=True,
                   method_standardize=None,
                   which_data="val",
                   first_id_plot=0, n_plot=None,
                   train_size=0.8, test_size=0.2, traintestrandomstate=42, traintestshuffle=True,
                   return_dictionary=False, **kwargs):
    if method_standardize is None:
        method_standardize = {"sfh": 0, "z": 0, "spectra": 0, "magnitudes": 0}
    else:
        for n in method_standardize.keys():
            if n not in ["sfh", "z", "spectra", "magnitudes"]:
                raise ValueError("method_standardize dictionary should only have the keys: "
                                 f"'sfh', 'z', 'spectra', 'magnitudes'. One of the keywords given is {n}")
    for key in ["sfh", "z", "spectra", "magnitudes"]:
        if key not in method_standardize.keys():
            method_standardize[key] = 0

    # Load Data
    print("[INFO] Loading data...")
    label_path = data_path.replace("Input_", "Label_")
    # metadata_path = data_path.replace("Input_", "MetaD_").replace(".fits", ".rda")

    if which_data not in ["all", "train", "val"]:
        raise ValueError(f"which_data needs to be 'all', 'train', or 'val'. given value: {which_data}")

    input_spectra, input_magnitudes, label_sfh, label_z, spectra_lambda, agevec, ageweight = \
        ERIK.loadfiles(input_path=data_path, labels_path=label_path,
                       method_standardize_label_sfh=method_standardize["sfh"],
                       method_standardize_label_z=method_standardize["z"],
                       method_standardize_spectra=method_standardize["spectra"],
                       method_standardize_magnitudes=method_standardize["spectra"])

    _, _, label_sfh_no_normalization, label_z_no_normalization, _, _, _ = \
        ERIK.loadfiles(input_path=data_path, labels_path=label_path,
                       method_standardize_label_sfh=0,
                       method_standardize_label_z=0,
                       method_standardize_spectra=0,
                       method_standardize_magnitudes=0)

    # Verified: If indices are added, they will STILL respect the output as if there were no indices
    split_train_test = train_test_split(input_spectra, input_magnitudes, label_sfh, label_z, label_sfh_no_normalization,
                                        label_z_no_normalization,
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

    saved_models = []
    for pqwe, model_p in enumerate(model_paths):
        saved_models.append(keras.models.load_model(model_p,
                                                    custom_objects={"smape_loss": XAVIER.Cerebro.smape_loss}))

    #################
    # Decide which data is going to be used
    if which_data == "all":
        # All the data
        in_data = np.concatenate([input_spectra, input_magnitudes], axis=1)
        idx = list(range(input_spectra.shape[0]))
    elif which_data == "train":
        # Only Training data
        in_data = np.concatenate([trainSpect, trainMag], axis=1)
        input_spectra = trainSpect
        input_magnitudes = trainMag
        label_sfh = trainLabSfh
        label_z = trainLabZ
        label_sfh_no_normalization = trainLabSfh_real
        label_z_no_normalization = trainLabZ_real
        idx = trainIndices
    elif which_data == "val":
        # Only Test data
        in_data = np.concatenate([testSpect, testMag], axis=1)
        input_spectra = testSpect
        input_magnitudes = testMag
        label_sfh = testLabSfh
        label_z = testLabZ
        label_sfh_no_normalization = testLabSfh_real
        label_z_no_normalization = testLabZ_real
        idx = testIndices
    else:
        raise ValueError

    if n_plot is None:
        n_plot = in_data.shape[0] - first_id_plot

    outputs = []
    for model in saved_models:
        outputs.append(model.predict(in_data[first_id_plot:(first_id_plot + n_plot)]))

    if return_loss:
        losses = []
        for model in saved_models:
            losses.append(
                [model.evaluate(np.array([in_data[idx, :]]),
                                [np.array([label_sfh[idx, :]]), np.array([label_z[idx, :]])],
                                batch_size=1, verbose=1)
                 for idx in range(first_id_plot, first_id_plot + n_plot)]
            )
    else:
        losses = None

    if return_dictionary:
        out = {"input_spectra": input_spectra, "input_magnitudes": input_magnitudes,
               "label_sfh": label_sfh, "label_z": label_z,
               "label_sfh_no_normalization": label_sfh_no_normalization,
               "label_z_no_normalization": label_z_no_normalization,
               "spectra_lambda": spectra_lambda, "agevec": agevec, "ageweight": ageweight,
               "saved_models": saved_models, "idx": idx, "outputs": outputs, "losses": losses}
        return out
    else:
        return input_spectra, input_magnitudes, label_sfh, label_z, label_sfh_no_normalization, \
               label_z_no_normalization, spectra_lambda, agevec, ageweight, saved_models, idx, outputs, losses


def model_colormap(model_paths, data_path, calculate_loss=False,
                   method_standardize=None,
                   which_data_is_going_to_be_used="val",
                   first_id_plot=0, n_plot=None, percentile=98,
                   mageburst=1e8, temp_folder="/Users/enrique/Documents/GitHub/LOGAN-SFH/tempFolder",
                   **kwargs):
    # Load Data
    if calculate_loss:
        return_loss = True
    else:
        return_loss = False
    loaded_data = evaluate_model(model_paths, data_path, return_loss=return_loss, method_standardize=method_standardize,
                                 which_data=which_data_is_going_to_be_used, first_id_plot=first_id_plot, n_plot=n_plot,
                                 **kwargs)

    (input_spectra, input_magnitudes, label_sfh, label_z, label_sfh_no_normalization, label_z_no_normalization,
     spectra_lambda, agevec, ageweight, saved_models, idx, outputs, losses) = loaded_data
    # remove tuple data
    del loaded_data

    # If n_plot is None, consider it as "all"
    if n_plot is None:
        n_plot = input_spectra.shape[0] - first_id_plot

    if calculate_loss:
        # Filter losses to keep only the relevant information
        combined_losses = [[k[0] for k in losses[i_model]] for i_model in range(len(saved_models))]
        log_combined_losses = [np.log10(_) for _ in combined_losses]

        metallicity_loss = [[k[1] for k in losses[i_model]] for i_model in range(len(saved_models))]
        log_metallicity_losses = [np.log10(_) for _ in metallicity_loss]

        sfh_losses = [[k[2] for k in losses[i_model]] for i_model in range(len(saved_models))]
        log_sfh_losses = [np.log10(_) for _ in sfh_losses]
    else:
        lick = []
        for imodel in range(len(model_paths)):
            lick.append(np.load(os.path.join(temp_folder, f"differences_{imodel}.npy")))
    # mass_weighted = label_sfh * ageweight
    # total_mass = np.sum(mass_weighted, axis=1)
    mass_weigted = label_sfh_no_normalization * ageweight
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

    if calculate_loss:
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
    else:
        for imodel in range(len(model_paths)):
            percentile_value = np.nanpercentile(lick[imodel],
                                             percentile, )
            id_percentile = np.logical_and(lick[imodel] < percentile_value, np.invert(np.isnan(lick[imodel])))
            # plot_colorbars(x=total_mass_no_norm[first_id_plot:(first_id_plot + n_plot)],
            #                y=mass_burst[first_id_plot:(first_id_plot + n_plot)],
            #                c=lick[imodel][first_id_plot:(first_id_plot + n_plot)],
            plot_colorbars(x=total_mass_no_norm[id_percentile],
                           y=mass_burst[id_percentile],
                           c=lick[imodel][id_percentile],
                           xlabel=r"$M_{total} (M_{\odot})$",
                           ylabel=r"$M_{burst} (M_{\odot}) (\leq0.1Gyr)$",
                           collabel=r'$Log_{10}(\ell)$',
                           title=r"Relative-MSE lick indexes based on $M_{total}$ and $M_{burst}$ for model #"
                                 + str(imodel),
                           savefig_name=f'/Users/enrique/Documents/GitHub/LOGAN-SFH/my_test_lick_{imodel}.png',
                           )
    print("------")


def verify_model(model_paths, data_path, name_models=None, return_loss=False, loss_function=None,
                 method_standardize=None,
                 which_data_is_going_to_be_used="val",
                 first_id_plot=0, n_plot=3,
                 **kwargs):
    # Load Data
    loaded_data = evaluate_model(model_paths, data_path, return_loss=return_loss,
                                 method_standardize=method_standardize,
                                 which_data=which_data_is_going_to_be_used, first_id_plot=first_id_plot, n_plot=n_plot,
                                 **kwargs)

    (input_spectra, input_magnitudes, label_sfh, label_z, label_sfh_real, label_z_real,
     spectra_lambda, agevec, ageweight, saved_models, idx, outputs, _) = loaded_data
    # remove tuple data
    del loaded_data

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
                                  f"{loss_function(label_sfh[i, :], k[0][(i - first_id_plot), :]).numpy():3f}")
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
                                  f"{loss_function(label_z[i, :], k[1][(i - first_id_plot), :]).numpy():4f}")
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
                   "ID": [idx[i]] * 5}
        tmp_id_names = {"ID": [idx[i]]}
        for q, k in enumerate(legend_name):
            if q == 0:
                continue
            else:
                tmp_id_names[f"name{q - 1}"] = [k]
        tmp_df = pd.DataFrame(tmp_dic)
        tmp_df.to_csv("/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file.pd", index=False)
        tmp_df2 = pd.DataFrame(tmp_sp)
        tmp_df3 = pd.DataFrame(tmp_mag)
        tmp_df4 = pd.DataFrame(tmp_id_names)
        tmp_df2.to_csv("/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file2.pd", index=False)
        tmp_df3.to_csv("/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file3.pd", index=False)
        tmp_df4.to_csv("/Users/enrique/Documents/GitHub/LOGAN-SFH/tmp_file4.pd", index=False)

        print(i)

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

        model_colormap(args.verify_model, data_path, model_names, **mainkwargs)
        # save_data_to_file(args.verify_model, data_path, **mainkwargs)
        # verify_model(args.verify_model, data_path, model_names, **arguments_standardize)

    if args.json_clean is not None:
        json_files = RAVEN.flatten_list(args.json_clean)
        for file in json_files:
            ERIK.prettyfy_json_file(file, args.verbose)


if __name__ == "__main__":
    # indices_lick = calculate_lick_indices_mass(calculate_for_input=True, verbosity_progress=400)
    # models_path = "/Users/enrique/Documents/GitHub/LOGAN-SFH/TrainedModels/MSE_reduced_wave_small_dataset/"
    models_path = "/Users/enrique/scpdata/2"
    datapath = os.path.join("/Users/enrique/Documents/GitHub/LOGAN-SFH/TrainingData/",
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

    model_colormap(models[:2], data_path=datapath, calculate_loss=False)


    # main({"verify_model": models[0:2], "data_path": datapath},
    #      model_names=filenames,
    #      method_standardize={"sfh": 5}, n_plot=10
    #      )





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
