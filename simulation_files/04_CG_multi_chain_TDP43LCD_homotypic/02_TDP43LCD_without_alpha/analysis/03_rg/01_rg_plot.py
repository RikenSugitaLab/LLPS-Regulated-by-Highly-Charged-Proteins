#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():

    n_temperature = 10
    n_frames = 2500

    data_out = open("TDP43_rg_vs_temperature.dat", "w")

    fig, axes = plt.subplots(1, 1, figsize=(6, 4))

    rg_ave_list = []
    rg_err_list = []
    for i in range(n_temperature):
        fname = "tdp43_rg_t{0:0>2d}_01.dat".format(i + 1)
        data_local = np.loadtxt(fname)

        rg_ave = np.mean(data_local)
        rg_err = np.std(data_local) / n_frames**0.5

        rg_ave_list.append(rg_ave)
        rg_err_list.append(rg_err)

        data_out.write("{0:>4d}    {1:>8.3f}    {2:>8.3f} \n".format(i * 10 + 260,   rg_ave, rg_err))

    X = [i * 10 + 260 for i in range(n_temperature)]
    axes.errorbar(X, rg_ave_list, yerr=rg_err_list, ls="-", c="b", elinewidth=1.5, ecolor="b", capsize=5, capthick=1.5, marker="o", mfc="b", mec="b", mew=1.5, ms=5)

    axes.set_xticks([20 * i + 260 for i in range(6)])
    axes.set_xticklabels([20 * i + 260 for i in range(6)], fontsize=10)
    axes.set_xlim(255, 355)

    axes.set_yticks([5 * i + 10 for i in range(10)])
    axes.set_yticklabels([5 * i + 10 for i in range(10)], fontsize=10)
    axes.set_ylim(24, 36)

    figname = "TDP43_rg_vs_temperature.svg"
    plt.savefig(figname, dpi=150)

if __name__ == '__main__':
    main()
