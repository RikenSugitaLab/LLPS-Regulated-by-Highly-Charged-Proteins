#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main(i_temperature):

    fig, axes = plt.subplots(1, 1, figsize=(9, 9))

    fname = "tdp43_contacts_t{0:0>2d}_01.dat".format(i_temperature)
    data_local = np.loadtxt(fname)
    zmax = np.max(data_local)
    print("Maximum of matrix:", np.max(data_local))

    axes.imshow(data_local, cmap="Reds", vmin=0.001, vmax=100, origin="lower")
    # axes.imshow(data_local, cmap="Blues", origin="lower")

    axes.set_xticks([0, 20, 40, 60, 80, 100, 120, 140, 160])
    axes.set_xticklabels([0, 20, 40, 60, 80, 100, 120, 140, 160], fontsize=10)
    axes.set_xlim(-1, 155)

    axes.set_yticks([0, 20, 40, 60, 80, 100, 120, 140, 160])
    axes.set_yticklabels([0, 20, 40, 60, 80, 100, 120, 140, 160], fontsize=10)
    axes.set_ylim(-1, 155)

    figname = "TDP43_T{0:3d}_intra_contact_map.svg".format(i_temperature * 10 + 260)
    plt.savefig(figname, dpi=150)

if __name__ == '__main__':
    for i in range(10):
        main(i + 1)
