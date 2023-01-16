#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():
    n_T = 1
    n_C = 9
    n_run = 5

    fname0 = "./tdp43_100_hero11_{0}_T_{1}_all_{2:0>2d}_density_averaged.dat"

    fig, ax = plt.subplots(1, 1, figsize=(9, 5), constrained_layout=True, sharex=False, sharey=False)
    # color_list = ["r", "g", "b", "m"]
    color_list = ["m"]
    X = [j * 10 + 10 for j in range(n_C)]

    for i_T in range(n_T):
        d_C_mean = []
        d_C_std  = []
        d_C_err  = []
        for i_C in range(n_C):
            d_data_all = []
            for i_run in range(n_run):
                fname = fname0.format(i_C * 10 + 10, i_T * 10 + 290, i_run + 1)
                d_data = np.loadtxt(fname, usecols=(4))
                d_data_all.append(d_data)
            d_data_flat = np.ravel(d_data_all)
            d_mean = d_data_flat.mean()
            d_std  = d_data_flat.std()
            d_C_mean.append(d_mean)
            d_C_std.append(d_std)
            d_C_err.append(d_std / d_data_flat.size**0.5)
            print(d_data_flat.size)

        # load homotypic density
        d0_data = np.loadtxt("/home/ctan/Port/HD1/HERO/TDP43/20220206_first_test/25_tdp43_100_multi_temperature_260_350_HPS-Urry_FIX_BONDLENGTH/analysis/density_evolution/tdp43_100_md2_T_{0:>02d}_density_averaged.dat".format(i_T + 4), usecols=(4))
        d0_mean = d0_data.mean() * 154

        # plot
        ax.errorbar(X, d_C_mean, yerr=d_C_err, ls="-", c=color_list[i_T], elinewidth=1.5, ecolor=color_list[i_T], capsize=5, capthick=1.5, marker="o", mfc="white", mec=color_list[i_T], mew=1.5, ms=5)
        ax.axhline(y=d0_mean, ls=":", lw=1, c=color_list[i_T])

        fout_fname = "tdp43_100_md2_T_{0:>02d}_particle_density.dat".format(i_T)
        fout = open(fout_fname, "w")
        fout.write("{0:>3d}    {1:>8.3f}    {2:>8.3f}   {3:>8.3f} \n".format(0, d0_mean, d0_data.std() * 154, 154 * d0_data.std() / d0_data.size ** 0.5 ))
        for iC in range(n_C):
            fout.write("{0:>3d}    {1:>8.3f}    {2:>8.3f}   {3:>8.3f} \n".format(iC * 10 + 10, d_C_mean[iC], d_C_std[iC], d_C_err[iC]))
        fout.close()

    # axis styles
    ax.set_xticks(X)
    ax.set_xticklabels(X, fontsize=12)
    ax.set_xlim(5, 95)
    ax.set_xlabel(r"$n_{Hero11}$", fontsize=16)
    ax.set_yticks([1.20, 1.30, 1.40, 1.50, 1.60])
    ax.set_yticklabels([1.20, 1.30, 1.40, 1.50, 1.60], fontsize=12)
    ax.set_ylim(1.11, 1.69)
    ax.set_ylabel(r"$\rho_{particle} (M)$", fontsize=16)

    figname = "TDP43_density_summary.svg".format()
    plt.savefig(figname, dpi=150)


if __name__ == '__main__':
    main()
