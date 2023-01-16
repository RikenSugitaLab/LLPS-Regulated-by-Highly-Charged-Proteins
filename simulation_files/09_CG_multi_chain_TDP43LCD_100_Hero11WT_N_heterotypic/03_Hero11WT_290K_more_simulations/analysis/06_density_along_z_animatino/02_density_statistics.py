#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():

    n_temperature = 4

    T = lambda i : i * 10 + 260

    n_pro1 = 100
    n_pro2 = [10 * i + 10 for i in range(9)]

    n_atom_pro1 = 154
    n_atom_pro2 = 99

    time_b = 3000
    time_e = 9000

    for i, n in enumerate(n_pro2):
        fig, axes = plt.subplots(n_temperature, 1, figsize=(9, 3 * n_temperature), constrained_layout=True, sharex=True, sharey=False)
        for i_t in range(n_temperature):
            temp = T(i_t)
            print("analyzing repo:", i, " : ", i_t)
            data1_fname = "../01_density_in_box/tdp43_100_hero11_{0}_T_{1}_tdp43_density_centered.dat".format(n, temp)
            data2_fname = "../01_density_in_box/tdp43_100_hero11_{0}_T_{1}_hero11_density_centered.dat".format(n, temp)
            pro1_density_data = np.loadtxt(data1_fname)
            pro2_density_data = np.loadtxt(data2_fname)

            pro1_density_histogram = np.mean(pro1_density_data[time_b:, :], axis=0)
            pro2_density_histogram = np.mean(pro2_density_data[time_b:, :], axis=0)

            X = [j for j in range(1, 101)]
            axes[i_t].plot(X, pro1_density_histogram, ls="-", c="b", lw=1.5)

            axes[i_t].set_xticks([10 * k + 1 for k in range(11)])
            axes[i_t].set_xticklabels([3 * k - 15 for k in range(11)], fontsize=12)
            axes[i_t].set_xlim(0, 100)
            axes[i_t].set_xlabel(r"z ($\times 100\AA$)", fontsize=16)

            axes[i_t].set_yticks([0.010 * k for k in range(11)])
            axes[i_t].set_yticklabels([10 * k for k in range(11)], fontsize=12)
            axes[i_t].set_ylim(0, 0.045)
            axes[i_t].set_ylabel(r"density of TDP43 (M)", fontsize=16)
            axes[i_t].tick_params(axis='y', labelcolor="b")

            ax_tmp = axes[i_t].twinx()
            ax_tmp.plot(X, pro2_density_histogram, ls="--", c="r", lw=1.5)

            ax_tmp.set_yticks([0.001 * k for k in range(11)])
            ax_tmp.set_yticklabels([1 * k for k in range(11)], fontsize=12)
            ax_tmp.set_ylim(0, 0.0085)
            ax_tmp.set_ylabel(r"density of Hero11 (M)", fontsize=16)
            ax_tmp.tick_params(axis='y', labelcolor="r")

        figname = "tdp43_hero11_{0}_averaged_density_along_z_t_{1}_{2}.svg".format(n, time_b, time_e)
        plt.savefig(figname, dpi=150)

if __name__ == '__main__':
    main()
