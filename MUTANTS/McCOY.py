# -*- coding: utf-8 -*-
# Measurements, COrrelation and sanitY of data

import argparse
import os
import sys


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


def verify_model(model_path, data_path, label_path):
    pass


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
        label_path = data_path.replace("Input_", "Label_")
        if args.verify_model[0] != "/":
            model_path = os.path.abspath(args.verify_model)
        else:
            model_path = args.verify_model

        # Verify Model
        verify_model(model_path, data_path, label_path)


if __name__ == "__main__":
    main({"verify_model": "model_file_best_GPUdeeper400.h5",
          "data_path": "20220119T154253_Hr3TXx"})
