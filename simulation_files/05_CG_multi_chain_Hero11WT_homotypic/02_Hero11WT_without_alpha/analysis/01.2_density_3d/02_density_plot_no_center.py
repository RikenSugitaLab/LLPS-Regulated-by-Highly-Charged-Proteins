#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():
    # ==========
    # parameters
    # ==========
    # trajectory name
    sim_name = "hero11_alpha_100_md1"

    # simulation temperatures
    T_list = [60 + i * 10 for i in range(30)]
    num_temperatures = len(T_list)

    # num of runs per temperature
    num_run = 5

    num_steps = 5000
    step_interval = 5
    plot_steps = num_steps // step_interval

    LOW_DENSITY_THRESHOLD = 0.003
    num_bins = 100

    # ======================================
    # density re-calculation and plotting...
    # ======================================
    fig, axes = plt.subplots(num_temperatures * num_run, 1, figsize=(9, num_temperatures * num_run), constrained_layout=True, sharex=True, sharey=False)
    for i, t in enumerate(T_list):
        print("Temperature:", t)
        for i_run in range(num_run):
            fname = "{0}_T_{1:0>2d}_{2:0>2d}.density.dat".format(sim_name, i + 1, i_run + 1)
            data_local = np.loadtxt(fname)[::step_interval, :]

            # ===========
            # plotting...
            # ===========
            axes[i * num_run + i_run].imshow(data_local.T, cmap="Reds", vmin=0.000001, vmax=.005, origin="lower")

            axes[i * num_run + i_run].set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
            axes[i * num_run + i_run].set_xlim(0, 1000)

            axes[i * num_run + i_run].set_yticks([0, 50, 100])
            axes[i * num_run + i_run].set_yticklabels([0, 10, 20], fontsize=10)
            axes[i * num_run + i_run].set_ylim(0, 100)
            axes[i * num_run + i_run].set_ylabel(r"z ($10^2\AA$)", fontsize=12)

            sub_leg = "T={0:>3d} K, run: {1:0>2d}".format(t, i_run + 1)
            axes[i * num_run + i_run].text(100, 75, sub_leg)

    axes[num_temperatures * num_run - 1].set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axes[num_temperatures * num_run - 1].set_xticklabels([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], fontsize=10)
    axes[num_temperatures * num_run - 1].set_xlim(0, 1000)
    axes[num_temperatures * num_run - 1].set_xlabel(r"MD steps ($\times 10^6$)", fontsize=12)

    figname = "{0}_T_density_no_center.png".format(sim_name)
    plt.savefig(figname, dpi=150)
    # plt.show()

if __name__ == '__main__':
    main()
