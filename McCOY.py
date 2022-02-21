# -*- coding: utf-8 -*-
# Measurements, COrrelation and sanitY of data

import argparse


def cleanlogfile(filename, filename_out = None):
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


if __name__ == "__main__":
    # Argument Parser (argparse)
    parser = argparse.ArgumentParser(description="Cleaning programm")

    parser.add_argument("-ln", "--logname", default=None, type=str,
                        help="If a filename is given, it will clean this log (i.e. remove extra characters")
    parser.add_argument("-clo", "--cleanlogoutput", default=None, type=str,
                        help="If a filename is given, will write the cleaned log in this location. "
                             "Requires a value in -ln to work.")
    args = parser.parse_args()
    if args.cleanlogoutput is not None and args.logname is None:
        parser.error("--cleanlogoutput requires --logname")

    if args.logname is not None:
        print(f"[INFO] Cleaning log '{args.logname}'")
        cleanlogfile(args.logname, args.cleanlogoutput)

    print("[LOG] FINISHED")
