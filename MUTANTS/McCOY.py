# -*- coding: utf-8 -*-
# Measurements, COrrelation and sanitY of data

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import keras

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
    print("min:", sum(lines2)/60)
    print("hours:", sum(lines2)/3600)

    lines.append(f"s_total: {sum(lines2)}")
    lines.append(f"m_total: {sum(lines2)/60}")
    lines.append(f"h_total: {sum(lines2)/3600}")

    lines = [line + "\n" for line in lines]
    with open(filename_out, "w+") as f:
        f.writelines(lines)


def verify_model(model_path, data_path):
    # Load Data
    print("[INFO] Loading data...")
    label_path = data_path.replace("Input_", "Label_")
    metadata_path = data_path.replace("Input_", "MetaD_").replace(".fits", ".json")
    input_spectra, input_magnitudes, label_sfh, label_z, spectra_lambda, agevec = \
        ERIK.loadfiles(input_path=data_path, labels_path=label_path)
    saved_model = keras.models.load_model(model_path, custom_objects={"smape_loss": XAVIER.Cerebro.smape_loss})
    in_data = np.concatenate([input_spectra, input_magnitudes], axis=1)

    # l_burst = []
    # for i in range(in_data.shape[0]):
    #     metadata = ERIK.getparametersfromid(metadata_path, i, verbose=0, returnfunction=False)
    #     if "mburst" in metadata.keys():
    #         l_burst.append(i)
    # print(l_burst)
    primer_valor = 30065
    output = saved_model.predict(in_data[primer_valor:(primer_valor + 1000)])
    for i in range(primer_valor, primer_valor + 1000):
        if label_sfh[i, 0] > 0.5:
            param = ERIK.getparametersfromid(metadata_path, i, verbose=0, returnfunction=False)
            if param["Zfinal"] == 0.02:
                print(param)
                print("mayor que 0.5")
                fig3, ax3 = plt.subplots()
                print("smape", XAVIER.Cerebro.smape_loss(output[0][i - primer_valor], label_sfh[i]))
                ax3.set_title(f"Test-SFH - {i} - {np.mean(XAVIER.Cerebro.smape_loss(output[0][i - primer_valor], label_sfh[i]))}")
                ax3.semilogx(agevec, label_sfh[i], "k-")
                ax3.semilogx(agevec, output[0][i - primer_valor], "g.-")
                plt.show()
                print("here")

    difference_sfh = np.subtract(output[0], label_sfh)
    difference_mag = np.subtract(output[1], label_z)

    fig1, ax1 = plt.subplots()
    ax1.set_title("Boxplot SFH")
    ax1.boxplot(difference_sfh)
    plt.show()
    fig2, ax2 = plt.subplots()
    ax2.set_title("Boxplot Z")
    ax2.boxplot(difference_mag)
    plt.show()

    for ids in range(10):

        metadata = ERIK.getparametersfromid("/Users/enrique/Documents/GitHub/LOGAN-SFH/MetadataOutput.json",
                                            ids, verbose=0, returnfunction=False)
        fig3, ax3 = plt.subplots()
        ax3.set_title(f"SFH - {ids}")
        ax3.semilogx(agevec, label_sfh[5000 + ids], "k-")
        ax3.semilogx(agevec, output[0][ids], "g.-")
        plt.savefig(f"/Users/enrique/Documents/GitHub/LOGAN-SFH/Test/SFH{ids}.png")
        fig4, ax4 = plt.subplots()
        ax4.set_title(f"Z - {ids}")
        ax4.semilogx(agevec, label_z[5000 + ids], "k-")
        ax4.semilogx(agevec, output[1][ids], "g.-")
        plt.savefig(f"/Users/enrique/Documents/GitHub/LOGAN-SFH/Test/Z{ids}.png")

    plt.show()
    print("test1234")


def main(update_values=None):
    # Argument Parser (argparse)
    parser = argparse.ArgumentParser(description="Cleaning programm")

    parser.add_argument("-ln", "--logname", default=None, type=str,
                        help="If a filename is given, it will clean this log (i.e. remove extra characters")
    parser.add_argument("-clo", "--cleanlogoutput", default=None, type=str,
                        help="If a filename is given, will write the cleaned log in this location. "
                             "Requires a value in --logname to work.")
    parser.add_argument("-vm", "--verify_model", default=None, type=str,
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
        print(f"[INFO] Verifying Model {args.verify_model} with Data:{args.data_path}")
        data_path = args.data_path
        if data_path[-5:] != ".fits":
            data_path = "Input_" + data_path + ".fits"
        if data_path[0] != "/":
            data_path = os.path.abspath(data_path)
        if args.verify_model[0] != "/":
            model_path = os.path.abspath(args.verify_model)
        else:
            model_path = args.verify_model

        # Verify Model
        verify_model(model_path, data_path)

    if args.json_clean is not None:
        json_files = RAVEN.flatten_list(args.json_clean)
        for file in json_files:
            ERIK.prettyfy_json_file(file, args.verbose)


if __name__ == "__main__":
    main({"verify_model": "/Users/enrique/Documents/GitHub/LOGAN-SFH/model_file_best_GPUdeeper400.h5",
          "data_path": "/Volumes/Elements/Outputs/Input_combined.fits"})
    #"data_path": "/Users/enrique/Documents/GitHub/LOGAN-SFH/OutputFolder/Input_combined.fits"})
    # main({"json_clean": [["/Users/enrique/Documents/GitHub/LOGAN-SFH/OutputFolder/Metadata_20220209T142129_fJnn24json",
    #                       "/Users/enrique/Documents/GitHub/LOGAN-SFH/OutputFolder/Metadata_20220210T122420_N9HRfM.json",
    #                       ]]
    #       })
