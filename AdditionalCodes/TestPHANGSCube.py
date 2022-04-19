from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def main():
    filename = "/Users/enrique/Documents/GitHub/LOGAN-SFH/TestPHANGS/ADP.2021-07-16T10_20_56.543.fits"
    filename = "/Users/enrique/Downloads/IC5332_MAPS_native.fits"

    hdul = fits.open(filename)
    # print(hdul)
    # hdul.info()
    hdr = hdul[0].header
    tmp = str(hdr)
    for i in range(len(hdr) + 1):
        print(tmp[(i*80):(i*80+80)])

    for i in range(1, len(hdul)):
        data = hdul[i].data
        data_ = data.copy()[~np.isnan(data)]
        vmin = np.percentile(data_, 5)
        vmax = np.percentile(data_, 95)
        fig, ax = plt.subplots()
        ax.imshow(data, vmin=vmin, vmax=vmax)
        ax.set_title(hdul[i].header["EXTNAME"])
        plt.savefig("/Users/enrique/Downloads/IC5332/" + str(i) + "_" + hdul[i].header["EXTNAME"] + ".png")
        print(hdul[i].header["EXTNAME"] + "-" + str(vmin) + "-" + str(vmax))
        pass
    print("Finished1")
    return

    for element in hdr.cards:
        print(element[0], "=", element[1], "///", element[2])
    print(bcolors.HEADER + "END" + bcolors.ENDC)

    tmp = str(hdr)
    for i in range(len(hdr) + 1):
        print(tmp[(i*80):(i*80+80)])
    print("Finished")









    hdul.close()


if __name__ == "__main__":
    main()
