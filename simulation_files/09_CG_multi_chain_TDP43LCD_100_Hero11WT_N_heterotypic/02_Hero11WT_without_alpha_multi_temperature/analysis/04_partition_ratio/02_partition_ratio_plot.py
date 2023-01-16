#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main(n_hero):
    temperature_list = [260 + i * 10 for i in range(4)]

    data_out_fname = "hero11_{0}_fraction_temperature.dat".format(n_hero)
    data_out = open(data_out_fname, "w")
    frac_ave_list = []
    frac_std_list = []
    for i, t in enumerate(temperature_list):
        data_fname = "./tdp43_100_hero11_{0}_t{1:>3d}.hero11.partition.dat".format(n_hero, t)
        partio_data = np.loadtxt(data_fname)
        n0 = partio_data.shape[0]
        X  = partio_data[n0 // 3:, 0]
        Y1 = partio_data[n0 // 3:, 1]
        Y2 = partio_data[n0 // 3:, 2]
        Y_sum = Y1 + Y2
        y_fraction = Y1 / Y_sum
        y_fraction_ave = np.mean(y_fraction)
        y_fraction_std = np.std(y_fraction)
        frac_ave_list.append(y_fraction_ave)
        # frac_std_list.append(y_fraction_std / len(Y1)**0.5)
        frac_std_list.append(y_fraction_std)
        data_out.write("{0:>4d}     {1:>8.3f}    {2:>8.3f}    {3:>8.3f} \n".format(t, y_fraction_ave, y_fraction_std, y_fraction_std / len(Y1)**0.5))

    fig, ax = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True, sharex=False, sharey=False)
    ax.errorbar(temperature_list, frac_ave_list, yerr=frac_std_list, ls="-", c="r", elinewidth=1.5, ecolor="r", capsize=5, capthick=1.5, marker="o", mfc="none", mec="r", mew=1.5, ms=5)

    ax.set_yticks([0.2 * i for i in range(6)])
    ax.set_yticklabels([0, "0.2", "0.4", "0.6", "0.8", "1.0"], fontsize=12)
    ax.set_ylim(0, 1)
    ax.set_ylabel("f(Hero11)", fontsize=16)

    plt.savefig("HERO_{0}_Partition_Ratio_vs_Temperature.svg".format(n_hero))

if __name__ == '__main__':
    for j in range(9):
        main(j * 10 + 10)
