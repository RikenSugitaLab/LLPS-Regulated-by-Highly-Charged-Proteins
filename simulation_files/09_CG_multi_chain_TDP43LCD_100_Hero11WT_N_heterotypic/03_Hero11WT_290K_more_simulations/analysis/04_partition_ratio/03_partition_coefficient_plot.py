#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():
    n_conc = 9
    temperature_list = [260 + i * 10 for i in range(4)]

    pc_ave_all = []
    pc_std_all = []
    pc_err_all = []
    for i, t in enumerate(temperature_list):
        data_out_fname = "hero11_T{0}_partition_coefficient.dat".format(t)
        data_out = open(data_out_fname, "w")
        pc_ave_list = []
        pc_std_list = []
        pc_err_list = []

        for j in range(n_conc):
            num_hero = (j + 1) * 10
            data_fname = "./tdp43_100_hero11_{0}_t{1:>3d}.hero11.partition.dat".format(num_hero, t)
            partio_data = np.loadtxt(data_fname)
            n0 = partio_data.shape[0]
            X = partio_data[n0 // 3:, 0]
            Y1 = partio_data[n0 // 3:, 1]
            Y2 = partio_data[n0 // 3:, 2]
            y_partition_coef = Y1 / Y2
            y_Pc_ave = np.mean(y_partition_coef)
            y_Pc_std = np.std(y_partition_coef)
            pc_ave_list.append(y_Pc_ave)
            pc_std_list.append(y_Pc_std)
            pc_err_list.append(y_Pc_std / len(Y1)**0.5)
            data_out.write("{0:>4d}     {1:>8.3f}    {2:>8.3f}    {3:>8.3f} \n".format(t, y_Pc_ave, y_Pc_std, y_Pc_std / len(Y1)**0.5))
        pc_ave_all.append(pc_ave_list[:])
        pc_std_all.append(pc_std_list[:])
        pc_err_all.append(pc_err_list[:])


    color_list = ["r", "g", "b", "m"]
    fig, ax = plt.subplots(1, 1, figsize=(6, 6), constrained_layout=True, sharex=False, sharey=False)
    for i, t in enumerate(temperature_list):
        X = [1*(j + 1) for j in range(n_conc)]
        Y = pc_ave_all[i]
        YERR = pc_std_all[i]
        ax.errorbar(X, Y, yerr=YERR, ls="-", c=color_list[i], elinewidth=1.5, ecolor=color_list[i], capsize=5, capthick=1.5, marker="o", mfc="none", mec=color_list[i], mew=1.5, ms=5)

    ax.set_yticks([0.5 * i for i in range(7)])
    ax.set_yticklabels([0, "0.5", "1.0", "1.5", "2.0", "2.5", "3.0"], fontsize=12)
    ax.set_ylim(0, 2.95)
    ax.set_ylabel("P", fontsize=16)
    # ax.set_yscale("log")

    ax.set_xlabel(r"$n_{Hero}$", fontsize=16)

    plt.savefig("HERO_Partition_coefficient.svg")

if __name__ == '__main__':
    main()
