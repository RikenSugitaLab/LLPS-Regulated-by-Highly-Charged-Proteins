#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main(i_temperature):

    fname = "tdp43_contacts_t{0:0>2d}_01.chain_pairwise.dat".format(i_temperature)
    data_local = np.loadtxt(fname)
    zmax = np.max(data_local)
    print("Maximum of matrix:", np.max(data_local))

    contact_num_1d = np.sum(data_local, axis=0)
    fout_name = "tdp43_contact_count_1d_T{0:>02d}.chain_pairwise.dat".format(i_temperature)
    fout = open(fout_name, "w")
    for j in range(len(contact_num_1d)):
        fout.write("{0:>3d}    {1:12.6f} \n".format(j + 1, contact_num_1d[j]))

    ###########################################################################
    #                                 plotting                                #
    ###########################################################################
    fig = plt.figure(figsize=(8, 8))

    # ================
    # make the 2d axis
    # ================
    # add the 2D contact map axis
    rect_2d = [0.1, 0.1, 0.6, 0.6]
    ax_2d = plt.axes(rect_2d, facecolor='w')
    ax_2d.set_xticks([20 * i - 1 for i in range(8)])
    ax_2d.set_xticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    ax_2d.set_xlim(-0.5, 153.5)
    ax_2d.set_xlabel("residue index", fontsize=16)
    ax_2d.set_xticks([20 * i - 1 for i in range(8)])
    ax_2d.set_xticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    ax_2d.set_xlim(-0.5, 153.5)
    ax_2d.set_xlabel("residue index", fontsize=16)

    # plot 2d contact probabilities
    ax_2d.imshow(data_local, cmap="Blues", vmin=0.001, vmax=0.4, origin="lower", aspect="equal")

    # ================
    # make the 1d axis
    # ================
    # add the 1d contact number count axis
    rect_1d = [0.1, 0.72, 0.6, 0.15]
    ax_1d = plt.axes(rect_1d, facecolor='w', sharex = ax_2d)

    # plot 1d contact probabilities
    X = [i for i in range(len(contact_num_1d))]
    Y = contact_num_1d
    # ax_1d.bar(X, Y, color="b")
    ax_1d.plot(X, Y, color="b")

    ax_1d.tick_params(axis="x", labelbottom=False)

    # plt.show()
    figname = "TDP43_T{0:3d}_inter_contact_map_and_contact_count.svg".format(i_temperature * 10 + 250)
    plt.savefig(figname, dpi=150)

if __name__ == '__main__':
    for i in range(10):
        main(i + 1)
