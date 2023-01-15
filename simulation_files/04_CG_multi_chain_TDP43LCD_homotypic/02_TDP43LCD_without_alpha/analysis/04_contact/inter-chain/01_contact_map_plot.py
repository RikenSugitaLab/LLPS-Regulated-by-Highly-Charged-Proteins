#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main(i_temperature):

    fig, axes = plt.subplots(1, 1, figsize=(9, 9))

    fname = "tdp43_contacts_t{0:0>2d}_01.dat".format(i_temperature)
    data_local = np.loadtxt(fname)
    zmax = np.max(data_local)
    print("Maximum of matrix:", np.max(data_local))

    axes.imshow(data_local, cmap="Blues", vmin=0.001, vmax=0.2, origin="lower")

    axes.set_xticks([20 * i - 1 for i in range(8)])
    axes.set_xticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    axes.set_xlim(-0.5, 153.5)

    axes.set_yticks([20 * i - 1 for i in range(8)])
    axes.set_yticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    axes.set_ylim(-0.5, 153.5)

    figname = "TDP43_T{0:3d}_inter_contact_map.svg".format(i_temperature * 10 + 250)
    plt.savefig(figname, dpi=150)

if __name__ == '__main__':
    for i in range(10):
        main(i + 1)
