import matplotlib

import argparse
import os
import sys
import shutil
# import csv
# import pickle
import codecs
# from typing import Callable, Any
from pylick.pylick.analysis import Galaxy
from pylick.pylick.indices import IndexLibrary
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.pylab as pl
import keras
from tensorflow.keras.utils import plot_model
# import tensorflow
from sklearn.model_selection import train_test_split
from math import nan


# sys.path.append("/Users/enrique/Documents/GitHub/LOGAN-SFH/MUTANTS")
# import MUTANTS.ERIK
# import MUTANTS.XAVIER
# import MUTANTS.McCOY
# import MUTANTS.RAVEN


def plot_colormaps(x, y, c, s=3, cmap='plasma_r', log="xy", xlim=None, ylim=None,
                   xlabel=None, ylabel=None, collabel=None, title=None, log_color=False,
                   savefig_name=None, figsize=(12, 9), same_axis=True, percentile_a=None, buffer=0.05, diagonal=False):
    fig, ax = plt.subplots(figsize=figsize)
    plt.title(title)
    if log_color:
        c = np.log10(c)
    sc = ax.scatter(x, y, c=c, s=s, cmap=cmap)
    if 'x' in log:
        ax.set_yscale('log')
    if "y" in log:
        ax.set_xscale('log')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label(collabel)
    if percentile_a is not None:
        percentile_min = np.nanpercentile(np.concatenate((x, y)), percentile_a * 100)
        percentile_max = np.nanpercentile(np.concatenate((x, y)), (1 - percentile_a) * 100)
        limits = (percentile_min - buffer * (percentile_max - percentile_min),
                  percentile_max + buffer * (percentile_max - percentile_min))
        xlim = limits
        ylim = limits
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    plt.tight_layout()
    if same_axis:
        x_limits = ax.get_xlim()
        y_limits = ax.get_ylim()
        limits = (min(x_limits[0], y_limits[0]), max(x_limits[1], y_limits[1]))
        ax.set_xlim(limits)
        ax.set_ylim(limits)
    if diagonal:
        x_tmp = ax.get_xlim()
        y_tmp = ax.get_ylim()
        ax.plot(x_tmp, y_tmp, "k-")
    if savefig_name is not None:
        print(f"[INFO] Saving image to: {savefig_name} ...")
        plt.savefig(savefig_name, facecolor='white', transparent=False, frameon=True)


def plot_lick_indices_comparison(in_spect, pre_spect, lick_indices_used,
                                 age_vec, lab_sfh, lab_z, pre_sfh, pre_z,
                                 wave_plot=np.arange(4700, 7501.25, 1.25),
                                 output_name="Lick_lines.png",
                                 force_error_sn=100, color_max=0.9, color_min=0.5,
                                 subplot_title_size=20, cmap=pl.cm.plasma_r,
                                 figsize=(60, 40), spectra_linewidth=4):
    lib = IndexLibrary()
    tex_names = lib.tex_names[lick_indices_used]
    n_lick = len(lick_indices_used)
    colors_map = cmap(np.linspace(0, 1, n_lick))
    fig = plt.figure(constrained_layout=True, figsize=figsize)
    gs = GridSpec(n_lick, 3, figure=fig)
    ax = []
    for il in range(n_lick):
        ax.append(fig.add_subplot(gs[il, 0]))
    # Spectra:
    ax.append(fig.add_subplot(gs[:-2 * int(n_lick / 4), 1:]))
    # SFH
    ax.append(fig.add_subplot(gs[-2 * int(n_lick / 4):-int(n_lick / 4), 1:]))
    # Z
    ax.append(fig.add_subplot(gs[-int(n_lick / 4):, 1:]))

    ax2 = []
    tmp1 = Galaxy("input", lick_indices_used, spec_wave=wave_plot,
                  spec_flux=in_spect,
                  spec_err=in_spect / force_error_sn,
                  spec_mask=None, meas_method='int',
                  plot=False, )
    tmp2 = Galaxy("predicted", lick_indices_used, spec_wave=wave_plot, spec_flux=pre_spect,
                  spec_err=pre_spect / force_error_sn,
                  spec_mask=None, meas_method='int',
                  plot=False, )

    for i_lick, lick in enumerate(lick_indices_used):
        regions = lib.regions[lick, :, :]
        xlimits = [regions[0, 0], regions[2, 1]]
        lick_filter = np.where(np.logical_and(xlimits[0] <= wave_plot, wave_plot <= xlimits[1]))
        ax[i_lick].set_title(
            fr"${tex_names[i_lick]}$ - in:{tmp1.vals[i_lick]:.3f} - pred:{tmp2.vals[i_lick]:.3f}",
            fontsize=subplot_title_size)
        ax[i_lick].plot(wave_plot[lick_filter], in_spect[lick_filter], color="black")
        if i_lick == n_lick - 1:
            ax[i_lick].set_xlabel(r"Wavelength $(\AA)$")
        ax[i_lick].set_ylabel("Input flux", color="black")
        ylimits = ax[i_lick].get_ylim()
        # Draw vertical lines
        for i1, color in enumerate(["blue", "black", "red"]):
            for i2 in range(2):
                ax[i_lick].axvline(x=regions[i1, i2], color=color)
            ax[i_lick].fill_between(wave_plot[lick_filter], 0, 1,
                                    where=np.logical_and(wave_plot[lick_filter] >= regions[i1, 0],
                                                         wave_plot[lick_filter] <= regions[i1, 1]),
                                    color=color, alpha=0.4)
        ax[i_lick].set_ylim(ylimits)
        # Right
        # append to ax2 list and read last element (avoids problems regarding reusing the same variable
        ax2.append(ax[i_lick].twinx())
        ax2[-1].plot(wave_plot[lick_filter], pre_spect[lick_filter], color="red")
        ax2[-1].set_ylabel("Predicted flux", color="red")

    # Spectra
    ax[-3].plot(wave_plot, in_spect, color="black", linewidth=spectra_linewidth)
    ax[-3].set_xlabel(r"Wavelength $(\AA)$")
    ax[-3].set_ylabel("Input flux", color="black")
    ax[-3].set_title("Spectra Comparison", fontsize=subplot_title_size)
    ylimits_spectra = ax[-3].get_ylim()
    # Draw vertical lines
    # Get specific values so that they do not overlap
    # (igual que las graficas en las que se comparan espectros y se normalizan a valores enteros)
    c_min = (ylimits_spectra[1] - ylimits_spectra[0]) * color_min + ylimits_spectra[0]
    c_max = (ylimits_spectra[1] - ylimits_spectra[0]) * color_max + ylimits_spectra[0]
    c_dis = (c_max - c_min) / n_lick
    for il, lick_ in enumerate(lick_indices_used):
        regions = lib.regions[lick_, :, :]
        lick_filter = np.where(np.logical_and(regions[0, 0] <= wave_plot, wave_plot <= regions[2, 1]))
        ax[-3].fill_between(wave_plot[lick_filter], c_max - il * c_dis, c_max - (il + 1) * c_dis,
                            color=colors_map[il], alpha=0.4)
        ax[-3].text((regions[0, 0] + regions[2, 1]) / 2, c_max - (il + 0.5) * c_dis, fr"${tex_names[il]}$",
                    color="black", horizontalalignment='center', verticalalignment='center')
    ax[-3].set_ylim(ylimits_spectra)
    ax3 = ax[-3].twinx()
    ax3.plot(wave_plot, pre_spect, color="red", linewidth=spectra_linewidth)
    ax3.set_ylabel("Predicted flux", color="red")

    # SFH
    ax[-2].plot(age_vec, lab_sfh, color="black", linewidth=spectra_linewidth)
    ax[-2].set_xscale('log')
    ax[-2].set_xlabel("LBT")
    ax[-2].set_ylabel("SFR", color="black")
    ax[-2].set_title("SFH Comparison", fontsize=subplot_title_size)
    ax2.append(ax[-2].twinx())
    ax2[-1].plot(age_vec, pre_sfh, color="red", linewidth=spectra_linewidth)
    ax2[-1].set_ylabel("SFR", color="red")

    # Z
    ax[-1].plot(age_vec, lab_z, color="black", linewidth=spectra_linewidth)
    ax[-1].set_xscale('log')
    ax[-1].set_xlabel("LBT")
    ax[-1].set_ylabel("Z", color="black")
    ax[-1].set_title("Z Comparison", fontsize=subplot_title_size)
    ax2.append(ax[-1].twinx())
    ax2[-1].plot(age_vec, pre_z, color="red", linewidth=spectra_linewidth)
    ax2[-1].set_ylabel("Z", color="red")

    plt.savefig(output_name, facecolor='white', transparent=False, frameon=True)


def main(imodel,
         models_path="/Users/enrique/Documents/GitHub/LOGAN-SFH/TrainedModels/MSLE_2/",
         temp_folder="/Users/enrique/Documents/GitHub/LOGAN-SFH/tempFolder_MSLE_2/",
         dummy_folder="/Users/enrique/Documents/GitHub/LOGAN-SFH/KK2/",
         n_spectra_to_plot=10,
         number_of_points_in_spectra_regenerated=1000,
         force_error_sn=100,
         percentile_a=0.01,
         ):

    if not os.path.exists(dummy_folder):
        os.mkdir(dummy_folder)

    lib = IndexLibrary()
    index_filter = list(range(number_of_points_in_spectra_regenerated))
    norm_mode = ["no_norm", "median", "mean"][0]
    datapath = os.path.join("/Users/enrique/Documents/GitHub/LOGAN-SFH/TrainingData/",
                            "Input_combined_wave_steps.fits")
    files = os.listdir(models_path)
    keywords_to_remove = ["epoch001", "imagen", "config"]
    for keyword in keywords_to_remove:
        files = [_ for _ in files if keyword not in _]
    files.sort()
    length_prefix_epoch = files[0].find("epoch")
    length_prefix_loss = files[0].find("loss")
    filenames = [f"e{_[(length_prefix_epoch + 5):(length_prefix_epoch + 8)]}-"
                 f"l{_[(length_prefix_loss + 4):-3]}" for _ in files]
    models = [os.path.join(models_path, f) for f in files if f[0] != "."]

    print("MODELS")
    print(models)

    loaded_model = keras.models.load_model(models[0])
    loaded_model.summary()
    plot_model(loaded_model, to_file=os.path.join(dummy_folder, "model_graph.png"),
               show_shapes=True, show_layer_names=True)
    del loaded_model

    name_input_spectra = os.path.join(temp_folder, "input_spectra.npy")
    name_input_spectra_steps = os.path.join(temp_folder, "input_spectra_steps.npy")
    name_input_spectra_regenerated = os.path.join(temp_folder, "input_spectra_regenerated.npy")
    name_wave = os.path.join(temp_folder, "wave.npy")
    name_files_predicted = [os.path.join(temp_folder, f"predict_spectra_{_}.npy") for _ in range(len(models))]
    name_files_predicted_step = [os.path.join(temp_folder, f"predict_spectra_{_}_step.npy") for _ in range(len(models))]
    out = {}

    input_spectra = np.load(name_input_spectra)
    input_spectra_steps = np.load(name_input_spectra_steps)
    input_spectra_regenerated = np.load(name_input_spectra_regenerated)
    input_spectra_number = input_spectra.shape[0]
    calculate_for_input = False
    verbosity_progress = 200
    # wave = np.load(name_wave)
    wave = np.arange(4700, 7501.25, 1.25)
    wave_step = (wave[:-1] + wave[1:]) / 2

    predicted_spectra = np.load(name_files_predicted[imodel])
    predicted_spectra_step = np.load(name_files_predicted_step[imodel])
    agevec = np.load(os.path.join(temp_folder, "agevec.npy"))
    ageweight = np.load(os.path.join(temp_folder, "ageweight.npy"))
    label_sfh_no_norm = np.load(os.path.join(temp_folder, "label_sfh_no_normalization.npy"))
    label_z_no_norm = np.load(os.path.join(temp_folder, "label_z_no_normalization.npy"))
    label_sfh = np.load(os.path.join(temp_folder, "label_sfh.npy"))
    label_z = np.load(os.path.join(temp_folder, "label_z.npy"))

    predicted_sfh = np.load(os.path.join(temp_folder, f"predict_sfh_{imodel}.npy"))
    predicted_z = np.load(os.path.join(temp_folder, f"predict_z_{imodel}.npy"))

    i_to_test = list(range(n_spectra_to_plot))
    xlimits = [4700, 7500]
    cut_waves = True
    buffer = 0.05

    os.makedirs(f"{dummy_folder}/ComparisonOfSpectra/model_{imodel}/", exist_ok=True)

    if input_spectra.shape[1] == 8200:
        input_spectra = input_spectra.T

    for i in i_to_test:
        print(f"[INFO] Plotting spectra for i: {i:3d}")
        fig, ax = plt.subplots(nrows=3, figsize=(12, 9))
        # Top
        # Left
        if cut_waves:
            wave_between_limits = np.where(np.logical_and(xlimits[0] <= wave, wave <= xlimits[1]))
            wave_cut = wave[wave_between_limits]
            y_max_l = np.max(input_spectra[i, wave_between_limits])
            y_min_l = np.min(input_spectra[i, wave_between_limits])
            xlim = (xlimits[0] - buffer * (xlimits[1] - xlimits[0]), xlimits[1] + buffer * (xlimits[1] - xlimits[0]))
            ylim_l = (y_min_l - buffer * (y_max_l - y_min_l), y_max_l + buffer * (y_max_l - y_min_l))
            y_max_r = np.max(predicted_spectra[i, wave_between_limits])
            y_min_r = np.min(predicted_spectra[i, wave_between_limits])
            ylim_r = (y_min_r - buffer * (y_max_r - y_min_r), y_max_r + buffer * (y_max_r - y_min_r))

        else:
            xlim = None
            ylim_l = None
            ylim_r = None
        ax[0].plot(wave_cut, input_spectra[i, :], color="black")
        ax[0].set_xlim(xlim)
        ax[0].set_ylim(ylim_l)
        ax[0].set_xlabel(r"Wavelength $(\AA)$")
        ax[0].set_ylabel("Input flux", color="black")
        # Right
        ax2 = ax[0].twinx()
        ax2.plot(wave, predicted_spectra[i, :], color="red")
        ax2.set_ylabel("Predicted flux", color="red")
        ax2.set_ylim(ylim_r)
        ############
        ax[1].plot(wave_step, input_spectra_steps[i, :], color="black")
        ax[1].set_xlim(xlim)
        # ax[1].set_ylim(ylim_l)
        ax[0].set_xlabel(r"Wavelength $(\AA)$")
        ax[1].set_ylabel("Input flux", color="black")
        # Right
        ax2 = ax[1].twinx()
        ax2.plot(wave_step, predicted_spectra_step[i, :], color="red")
        ax2.set_ylabel("Predicted flux", color="red")
        # ax2.set_ylim(ylim_r)
        # Bottom
        # Left
        ax[2].plot(agevec, label_sfh[i, :], color="black")
        ax[2].set_xscale('log')
        ax[2].set_xlabel("LBT")
        ax[2].set_ylabel("SFR", color="black")
        # Right
        ax4 = ax[2].twinx()
        ax4.plot(agevec, predicted_sfh[i, :], color="red")
        ax4.set_ylabel("SFR", color="red")
        plt.tight_layout()
        plt.savefig(f"{dummy_folder}/ComparisonOfSpectra/model_{imodel}/SpectraComparison_{i}.png",
                    facecolor='white', transparent=False, frameon=True)

    # index_list = np.array([53, 54, 55, 57, 58])         # Only interested lines
    index_list = np.arange(53, 66)  # All indices in range
    names = lib.names[index_list]
    tex_names = lib.tex_names[index_list]
    input_spectra_number = 5

    if not os.path.exists(os.path.join(dummy_folder, "lines")):
        os.mkdir(os.path.join(dummy_folder, "lines"))

    if not os.path.exists(os.path.join(dummy_folder, "lines", f"lick_{imodel}")):
        os.mkdir(os.path.join(dummy_folder, "lines", f"lick_{imodel}"))

    for i_spectra in range(n_spectra_to_plot):
        plot_lick_indices_comparison(in_spect=input_spectra_regenerated[i_spectra, :],
                                     pre_spect=predicted_spectra[i_spectra, :], lick_indices_used=index_list,
                                     age_vec=agevec, lab_sfh=label_sfh[i_spectra, :], lab_z=label_z[i_spectra, :],
                                     pre_sfh=predicted_sfh[i_spectra, :], pre_z=predicted_z[i_spectra, :],
                                     wave_plot=wave, force_error_sn=force_error_sn,
                                     output_name=os.path.join(dummy_folder, "lines", f"lick_{imodel}",
                                                              f"Lick_lines_{i_spectra}.png"))

    output_name_input_spectra_lick = "lick_idx_input_spectra_regenerated_"
    input_spectra = np.load(name_input_spectra_regenerated)
    input_spectra_number = min(input_spectra.shape[0], number_of_points_in_spectra_regenerated)
    index_values_spectra = np.zeros((input_spectra_number, len(index_list)), dtype=float)

    # Normalization
    if norm_mode == "no_norm":
        pass
    elif norm_mode == "mean":
        input_spectra = input_spectra / np.mean(input_spectra, axis=1)[:, None]
    elif norm_mode == "median":
        input_spectra = input_spectra / np.median(input_spectra, axis=1)[:, None]

    print("[INFO] Calculating lick indices for input spectra...")
    start_time = datetime.now()
    for i in range(input_spectra_number):
        if i % verbosity_progress == 0:
            print(f"Current case for input spectra : {i}/{input_spectra_number}")
        tmp = Galaxy("Input_" + str(i), index_list, spec_wave=wave, spec_flux=input_spectra[i, :],
                     spec_err=input_spectra[i, :] / force_error_sn, spec_mask=None, meas_method='int',
                     plot=False, plot_settings={'outpath': temp_folder + "/Graficas/"}, z=0.001)
        index_values_spectra[i, :] = tmp.vals
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    print("[INFO] Finished calculating lick indices for input spectra.")
    print(f"[INFO] storing in {os.path.join(temp_folder, output_name_input_spectra_lick + norm_mode + '.npy')} ...")
    np.save(os.path.join(temp_folder, output_name_input_spectra_lick + norm_mode + '.npy'), index_values_spectra)

    predicted_spectra_name = os.path.join(temp_folder, f"predict_spectra_{imodel}.npy")
    predicted_spectra = np.load(predicted_spectra_name)
    output_name_predicted_spectra_lick = "lick_idx_predicted_"
    predicted_spectra_number = min(predicted_spectra.shape[0], number_of_points_in_spectra_regenerated)
    index_values_predicted = np.zeros((predicted_spectra_number, len(index_list)), dtype=float)

    # Normalization
    if norm_mode == "no_norm":
        pass
    elif norm_mode == "mean":
        predicted_spectra = predicted_spectra / np.mean(predicted_spectra, axis=1)[:, None]
    elif norm_mode == "median":
        predicted_spectra = predicted_spectra / np.median(predicted_spectra, axis=1)[:, None]

    print(f"[INFO] Calculating lick indices for predicted spectra for model #{imodel}...")
    start_time = datetime.now()
    for i in range(predicted_spectra_number):
        if i % verbosity_progress == 0:
            print(f"Current case for predicted spectra model #{imodel} : {i:4d}/{predicted_spectra_number}")
        # FIXME: There are some cases where there are NaNs (maybe because sfh/z is predicted as <0?)
        if not np.isnan(predicted_spectra[i, 0]):
            tmp = Galaxy("Predicted_" + str(i), index_list, spec_wave=wave, spec_flux=predicted_spectra[i, :],
                         spec_err=predicted_spectra[i, :] / force_error_sn, z=0.001, spec_mask=None, meas_method='int')
            index_values_predicted[i, :] = tmp.vals
        else:
            index_values_predicted[i, :] = [nan] * len(index_list)
    end_time = datetime.now()
    print('Duration: {}'.format(end_time - start_time))
    print(f"[INFO] Finished calculating lick indices for predicted spectra for model #{imodel}.")
    print(
        f"[INFO] storing in {os.path.join(temp_folder, output_name_predicted_spectra_lick + norm_mode + '_' + str(imodel) + '.npy')} ...")
    np.save(os.path.join(temp_folder, output_name_predicted_spectra_lick + norm_mode + '_' + str(imodel) + '.npy'),
            index_values_predicted)

    mageburst = 1e8
    mass_weigted = label_sfh_no_norm * ageweight
    total_mass_no_norm = np.sum(mass_weigted, axis=1)
    mass_only_burst_idx = agevec < mageburst
    mass_burst = np.sum(mass_weigted[:, mass_only_burst_idx], axis=1)
    burst_ratio = np.log10(mass_burst / total_mass_no_norm)

    print(norm_mode)

    lick_spectra = np.load(os.path.join(temp_folder, f"lick_idx_input_spectra_regenerated_{norm_mode}.npy"))
    lick_predict = np.load(os.path.join(temp_folder, f"lick_idx_predicted_{norm_mode}_{imodel}.npy"))
    if not os.path.exists(os.path.join(dummy_folder, "lines", f"lick_{imodel}", "index_comparison")):
        os.mkdir(os.path.join(dummy_folder, "lines", f"lick_{imodel}", "index_comparison"))

    if not os.path.exists(os.path.join(dummy_folder, "lines", f"lick_{imodel}")):
        os.mkdir(os.path.join(dummy_folder, "lines", f"lick_{imodel}"))
    if not os.path.exists(os.path.join(dummy_folder, "lines", f"lick_{imodel}/index_comparison")):
        os.mkdir(os.path.join(dummy_folder, "lines", f"lick_{imodel}/index_comparison"))

    for iline, name, tex_name in zip(range(len(names)), names, tex_names):
        true_index = lick_spectra[:, iline]
        predicted_index = lick_predict[:, iline]
        plot_colormaps(x=true_index[index_filter],
                       y=predicted_index[index_filter],
                       c=burst_ratio[index_filter],
                       xlabel=fr"True value for ${tex_name}$",
                       ylabel=fr"Predicted value for ${tex_name}$",
                       collabel=r'$Log_{10}(M_{burst}/M_{total}) (Burst \equiv \sum M_{age\leq0.1Gyr})$',
                       title=fr"Lick indices comparison true vs predicted [${tex_name}$ line] - model: #{imodel}",
                       log_color=False,
                       log="",
                       percentile_a=percentile_a,
                       savefig_name=f'{dummy_folder}/lines/lick_{imodel}/index_comparison/'
                                    f'lick_{name}_{norm_mode}{f"_{percentile_a}" if percentile_a is not None else ""}.png',
                       diagonal=True
                       )


def sort_images(list_of_folders):
    for folder in list_of_folders:
        print(folder)
        subfolders = os.listdir(os.path.join(folder, "lines"))
        os.makedirs(os.path.join(folder, "line_evolution"), exist_ok=True)
        for i in range(10):
            os.makedirs(os.path.join(folder, "line_evolution", str(i)), exist_ok=True)
        for subfolder in subfolders:
            print(subfolder)
            if subfolder == ".DS_Store":
                continue
            for i in range(10):
                os.symlink(os.path.join(folder, "lines", subfolder,
                                        f"Lick_lines_{i}.png"),
                           os.path.join(folder, "line_evolution", str(i),
                                        f"{i}_{int(subfolder.split('_')[1]):02d}.png"))




if __name__ == "__main__":
    # models_folder_path = "/Users/enrique/Documents/GitHub/LOGAN-SFH/TrainedModels/MSLE_1/"
    # temp_folder_path = "/Users/enrique/Documents/GitHub/LOGAN-SFH/tempFolder_MSLE_1/"
    # dummy_folder_path = "/Users/enrique/Documents/GitHub/LOGAN-SFH/KK1/"
    # for imodel in range(7):
    #     main(imodel=imodel,
    #          models_path=models_folder_path,
    #          temp_folder=temp_folder_path,
    #          dummy_folder=dummy_folder_path,
    #          n_spectra_to_plot=10,
    #          number_of_points_in_spectra_regenerated=1000,
    #          force_error_sn=100,
    #          percentile_a=0.01,
    #          )
    sort_images(["/Users/enrique/Documents/GitHub/LOGAN-SFH/KK1/"])
    # sort_images(["/Users/enrique/Documents/GitHub/LOGAN-SFH/KK/",
    #              "/Users/enrique/Documents/GitHub/LOGAN-SFH/KK2/",
    #              "/Users/enrique/Documents/GitHub/LOGAN-SFH/KK3/",
    #              "/Users/enrique/Documents/GitHub/LOGAN-SFH/KK4/",
    #              "/Users/enrique/Documents/GitHub/LOGAN-SFH/KK5/"])
    # models_folder_path = "/Users/enrique/Documents/GitHub/LOGAN-SFH/TrainedModels/MSLE_5/"
    # temp_folder_path = "/Users/enrique/Documents/GitHub/LOGAN-SFH/tempFolder_MSLE_5/"
    # dummy_folder_path = "/Users/enrique/Documents/GitHub/LOGAN-SFH/KK5/"
    # for imodel in range(17):
    #     main(imodel=imodel,
    #          models_path=models_folder_path,
    #          temp_folder=temp_folder_path,
    #          dummy_folder=dummy_folder_path,
    #          n_spectra_to_plot=10,
    #          number_of_points_in_spectra_regenerated=1000,
    #          force_error_sn=100,
    #          percentile_a=0.01,
    #          )