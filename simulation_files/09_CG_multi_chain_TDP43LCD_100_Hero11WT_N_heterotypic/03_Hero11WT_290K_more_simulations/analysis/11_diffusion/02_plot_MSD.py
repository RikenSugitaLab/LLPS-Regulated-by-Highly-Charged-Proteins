#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main(T, num_hero, pro_name):
    n_run = 5

    data_fname0 = "./tdp43_100_hero11_{0:d}_T{1:>3d}_md3_r{2:0>2d}.{3}_MSD.dat"

    n_dt = 6
    dt_all = np.zeros((n_dt, n_run), dtype=int)
    dense_MSD_data_all = np.zeros((n_dt, n_run), dtype=int)
    dilute_MSD_data_all = np.zeros((n_dt, n_run), dtype=int)
    for j in range(n_run):
        MSD_data_local = np.loadtxt(data_fname0.format( num_hero, T, j + 1, pro_name))
        dt_all[:, j] = MSD_data_local[:, 0]
        dense_MSD_data_all[:, j] = MSD_data_local[:, 1]
        dilute_MSD_data_all[:, j] = MSD_data_local[:, 2]

    X = dt_all[:, 0]
    dense_MSD_ave = np.mean(dense_MSD_data_all, axis=1)
    dense_MSD_std = np.std(dense_MSD_data_all, axis=1)
    dilute_MSD_ave = np.mean(dilute_MSD_data_all, axis=1)
    dilute_MSD_std = np.std(dilute_MSD_data_all, axis=1)

    fig, ax = plt.subplots(2, 1, figsize=(9, 8), constrained_layout=True, sharex=True, sharey=False)
    ax[0].errorbar(X, dense_MSD_ave, yerr=dense_MSD_std, ls="-", c="r", elinewidth=1.5, ecolor="r", capsize=5, capthick=1.5, marker="o", mfc="white", mec="r", mew=1.5, ms=5)
    ax[1].errorbar(X, dilute_MSD_ave, yerr=dilute_MSD_std, ls="-", c="b", elinewidth=1.5, ecolor="b", capsize=5, capthick=1.5, marker="o", mfc="white", mec="b", mew=1.5, ms=5)

    figname = "tdp43_100_hero11_{1:d}_T{0:3d}_{2}_MSD.png".format(T, num_hero, pro_name)
    plt.savefig(figname, dpi=150)

    # ======================
    # fit diffusion constant
    # ======================
    dense_fit_coef,  res = np.polyfit(X, dense_MSD_ave,  1)
    dilute_fit_coef, res = np.polyfit(X, dilute_MSD_ave, 1)
    data_fname0 = "./tdp43_100_hero11_{1}_T_{0:>02d}.{2}_D.dat".format(T, num_hero, pro_name)
    data_of = open(data_fname0, "w")
    data_of.write("{0:18.12f}       {1:18.12f}  \n".format(dense_fit_coef / 6.0, dilute_fit_coef / 6.0))

if __name__ == '__main__':
    for j in range(9):
        main(290, j * 10 + 10, "tdp43")
        main(290, j * 10 + 10, "hero11")
