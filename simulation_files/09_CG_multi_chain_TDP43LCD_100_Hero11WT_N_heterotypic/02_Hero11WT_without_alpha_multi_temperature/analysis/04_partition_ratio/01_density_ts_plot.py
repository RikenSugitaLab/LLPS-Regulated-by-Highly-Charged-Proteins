#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


def main(num_hero):
    temp_list = [260 + i * 10 for i in range(4)]
    n_temperature = len(temp_list)

    system_name = "tdp43_100_hero11_{0}".format(num_hero)

    # plot...
    fig, axes = plt.subplots(n_temperature, 1, figsize=(9, 2*n_temperature), constrained_layout=True, sharex=True, sharey=False)
    for it, temperature in enumerate(temp_list):
        hero_data = np.loadtxt("{0}_t{1:>03d}.hero11.partition.dat".format(system_name, temperature))[::3, :]
        X = hero_data[:, 0]
        Y1 = hero_data[:, 1]
        Y2 = hero_data[:, 2]
        Y_sum = Y1 + Y2
        y1 = Y1 / Y_sum
        y2 = Y2 / Y_sum
        axes[it].plot(X, y1, c="r", lw=2)  # in
        axes[it].plot(X, y2, c="g", lw=2)  # out

        axes[it].set_yticks([k *0.2 for k in range(6)])
        axes[it].set_yticklabels(["0", "0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=12)
        axes[it].set_ylim(0, 1)
        axes[it].set_ylabel("fraction", fontsize=16)
        axes[it].set_xticks([k * 500 for k in range(7)])
        axes[it].set_xticklabels([k * 15 for k in range(7)], fontsize=12)
        axes[it].set_xlim(0, 3000)
        axes[it].set_xlabel(r"time steps ($\times 10^6$)", fontsize=16)
        axes[it].grid(axis="y", alpha=0.2)

    figname = "{0}_partition_timeseries.png".format(system_name)
    plt.savefig(figname, dpi=150)

if __name__ == '__main__':
    for j in range(9):
        main(10 * (j + 1))
