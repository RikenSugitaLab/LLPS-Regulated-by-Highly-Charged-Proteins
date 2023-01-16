#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():
    n_T = 4
    n_C = 9

    fT = lambda i : i * 10 + 260
    fC = lambda i : i * 10 + 10

    fig, axes = plt.subplots(n_C, n_T, figsize=(4 * n_T, 3 * n_C), constrained_layout=True, sharex=True, sharey=True)

    for i_T in range(n_T):
        for i_C in range(n_C):
            print("analyzing #T: ", i_T, "  #Hero: ", i_C)
            data_1_fname = "./md3/tdp43_100_hero11_{0}_md3_T_{1}.tdp43.density.dat".format(fC(i_C), fT(i_T))
            data_2_fname = "./md3/tdp43_100_hero11_{0}_md3_T_{1}.hero11.density.dat".format(fC(i_C), fT(i_T))

            # ==================
            # histogram of TDP43
            # ==================
            data_1 = np.loadtxt(data_1_fname)
            data_1_nonzero = data_1[np.nonzero(data_1)]
            d1_hist, d1_edge = np.histogram(data_1_nonzero, bins=100, density=True)
            X1 = 0.5 * (d1_edge[:-1] + d1_edge[1:])
            axes[i_C, i_T].plot(X1, d1_hist, "b-")

            # =================
            # histogram of Hero
            # =================
            data_2 = np.loadtxt(data_2_fname)
            data_2_nonzero = data_2[np.nonzero(data_2)]
            d2_hist, d2_edge = np.histogram(data_2_nonzero, bins=100, density=True)
            X2 = 0.5 * (d2_edge[:-1] + d2_edge[1:])
            axes[i_C, i_T].plot(X2, d2_hist, "r-")

            axes[i_C, i_T].set_yscale("log")

            if i_C == n_C - 1:
                axes[i_C, i_T].set_xlabel("Density (M)", fontsize=16)
            if i_T == 0:
                axes[i_C, i_T].set_ylabel("Histogram", fontsize=16)

    # plt.show()
    plt.savefig("Density_distributions.png", dpi=150)


if __name__ == '__main__':
    main()
