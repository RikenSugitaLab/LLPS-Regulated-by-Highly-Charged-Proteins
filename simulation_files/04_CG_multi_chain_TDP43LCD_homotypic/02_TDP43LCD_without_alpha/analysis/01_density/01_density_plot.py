#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():
    T_list = [260 + i * 10 for i in range(10)]
    num_temperatures = len(T_list)
    num_run = 1
    num_steps = 5000
    step_interval = 5

    fig, axes = plt.subplots(num_temperatures, 1, figsize=(9, num_temperatures), constrained_layout=True, sharex=True, sharey=False)
    for i, t in enumerate(T_list):
        print("Temperature:", t)
        # fname = "tdp43_md1_T_{0:0>2d}.density.dat".format(i + 1)
        fname = "tdp43_md2_T_{0:0>2d}.density.dat".format(i + 1)
        data_local = np.loadtxt(fname)[::step_interval, :]

        new_index = [i for i in range(50, 99)] + [i for i in range(50)]
        data_shift = data_local[:, new_index]

        axes[i].imshow(data_shift.T, cmap="Blues", vmin=0.000001, vmax=.005, origin="lower")

        axes[i].set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
        axes[i].set_xlim(0, 1000)

        axes[i].set_yticks([0, 50, 100])
        axes[i].set_yticklabels([0, 15, 30], fontsize=10)
        axes[i].set_ylim(0, 100)
        axes[i].set_ylabel(r"z ($10^2\AA$)", fontsize=12)

        sub_leg = "T={0:>3d} K".format(t)
        axes[i].text(100, 75, sub_leg)

    axes[num_temperatures - 1].set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axes[num_temperatures - 1].set_xticklabels([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], fontsize=10)
    axes[num_temperatures - 1].set_xlim(0, 1000)
    axes[num_temperatures - 1].set_xlabel(r"MD steps ($\times 10^6$)", fontsize=12)

    # figname = "tdp43_md1_T_density.png".format(t)
    figname = "tdp43_md2_T_density.png".format(t)
    plt.savefig(figname, dpi=150)

if __name__ == '__main__':
    main()
