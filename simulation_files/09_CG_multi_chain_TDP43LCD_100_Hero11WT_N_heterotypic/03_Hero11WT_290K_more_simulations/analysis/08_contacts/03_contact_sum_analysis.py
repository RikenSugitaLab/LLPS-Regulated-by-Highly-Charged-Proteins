#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


def main(phase_state):
    T_list = [290]
    n_list = [10 * i + 10 for i in range(9)]

    n_run = 5

    name_pro1 = "tdp43"
    name_pro2 = "hero11"

    fig, axes = plt.subplots(2, 1, figsize=(5, 5), constrained_layout=True, sharex=True, sharey=False)

    data_name = "./md3/tdp43_100_hero11_{0}_T{1}_md3_r{2:0>2d}.tdp43-{3}.{4}.contact_count_ts.dat"
    color_list = ["r", "g", "b", "m"]
    for i_T, T in enumerate(T_list):
        output_fname = "tdp43_hero11_{0}_contact_count_sum_T_{1:03d}.dat".format(phase_state, T)
        of_tmp = open(output_fname, "w")
        Y1_mean = []
        Y2_mean = []
        Ys_mean = []
        Y1_std = []
        Y2_std = []
        Ys_std = []
        Y1_err = []
        Y2_err = []
        Ys_err = []
        for n in n_list:
            y1 = []
            y2 = []
            for i_run in range(n_run):
                fname = data_name.format(n, T, i_run + 1, name_pro1, phase_state)
                with open(fname, "r") as fin:
                    for line in fin:
                        words = line.split()
                        if len(words) > 0:
                            for w in words:
                                if float(w) > -1e-6:
                                    y1.append(float(w))
                fname = data_name.format(n, T, i_run + 1, name_pro2, phase_state)
                with open(fname, "r") as fin:
                    for line in fin:
                        words = line.split()
                        if len(words) > 0:
                            for w in words:
                                if float(w) > -1e-6:
                                    y2.append(float(w))
            y1_m = np.mean(y1)
            y1_s = np.std(y1)
            y1_e = y1_s / len(y1)**0.5
            y2_m = np.mean(y2)
            y2_s = np.std(y2)
            y2_e = y2_s / len(y2)**0.5
            ys_m = y1_m + y2_m
            ys_s = (y1_s**2 + y2_s**2)**0.5
            ys_e = y1_e + y2_e
            Y1_mean.append(y1_m)
            Y2_mean.append(y2_m)
            Ys_mean.append(ys_m)
            Y1_std.append(y1_s)
            Y2_std.append(y2_s)
            Ys_std.append(ys_s)
            Y1_err.append(y1_e)
            Y2_err.append(y2_e)
            Ys_err.append(ys_e)

        axes[0].errorbar(n_list, Y1_mean, yerr=Y1_err, ls="--", c=color_list[i_T], elinewidth=1.5, ecolor=color_list[i_T], capsize=5, capthick=1.5, marker="o", mfc="white", mec=color_list[i_T], mew=1.5, ms=5, lw=1)
        axes[1].errorbar(n_list, Y2_mean, yerr=Y2_err, ls="-.", c=color_list[i_T], elinewidth=1.5, ecolor=color_list[i_T], capsize=5, capthick=1.5, marker="o", mfc="white", mec=color_list[i_T], mew=1.5, ms=5, lw=1)
        axes[0].errorbar(n_list, Ys_mean, yerr=Ys_err, ls="-", c=color_list[i_T], elinewidth=1.5, ecolor=color_list[i_T], capsize=5, capthick=1.5, marker="o", mfc="white", mec=color_list[i_T], mew=1.5, ms=5, lw=2)
        outstring = "{0:3d}    {1:12.4f} {2:12.4f} {3:12.4f}   {4:12.4f} {5:12.4f} {6:12.4f}   {7:12.4f} {8:12.4f} {9:12.4f} \n"
        for j in range(9):
            of_tmp.write(outstring.format(n_list[j], Y1_mean[j], Y1_std[j], Y1_err[j], Y2_mean[j], Y2_std[j], Y2_err[j], Ys_mean[j], Ys_std[j], Ys_err[j]))
        of_tmp.close()

    # add TDP43 homotypic LLPS data
    o_data_name = "/home/ctan/Port/HD1/HERO/TDP43/20220206_first_test/25_tdp43_100_multi_temperature_260_350_HPS-Urry_FIX_BONDLENGTH/analysis/contact_map/inter-chain/tdp43_100_t{0:0>3d}_md2.tdp43-tdp43.{1}.contact_count_ts.dat"
    output_fname = "tdp43_homo_{0}_contact_count_sum_T_all.dat".format(phase_state)
    of_tmp = open(output_fname, "w")
    for i_T, T in enumerate(T_list):
        fname = o_data_name.format(i_T + 4, phase_state)
        # cnt_11 = np.loadtxt(fname, usecols=(1))
        # cnt_sum_1 = np.sum(cnt_11)
        y_all = []
        with open(fname, "r") as fin:
            for line in fin:
                words = line.split()
                if len(words) > 0:
                    for w in words:
                        cn = float(w)
                        if cn >= 0:
                            y_all.append(cn)
        y_arr = np.array(y_all)
        y_mean = np.mean(y_arr)
        y_std  = np.std(y_arr)
        y_err  = y_std / len(y_all)**0.5
        axes[0].axhline(y=y_mean, ls=":", lw=1, c=color_list[i_T])
        axes[0].axhspan(y_mean - y_err, y_mean + y_err, color=color_list[i_T], alpha=0.3)

        outstring = "{0:3d}    {1:12.4f} {2:12.4f} {3:12.4f} \n"
        of_tmp.write(outstring.format(T, y_mean, y_std, y_err))

    # axis styles
    axes[0].set_xticks(n_list)
    # axes[0].set_xticklabels(n_list, fontsize=12)
    axes[0].set_xlim(5, 95)
    # axes[0].set_xlabel(r"$n_{Hero11}$", fontsize=16)
    # axes[0].set_yticks([2000, 3000, 4000, 5000])
    # axes[0].set_yticklabels([2, 3, 4, 5], fontsize=12)
    # axes[0].set_ylim(1950, 5050)
    axes[0].set_ylabel(r"$N_{contacts} (\times 1000)$", fontsize=16)

    axes[1].set_xticks(n_list)
    axes[1].set_xticklabels(n_list, fontsize=12)
    axes[1].set_xlim(5, 95)
    axes[1].set_xlabel(r"$n_{Hero11}$", fontsize=16)
    # axes[1].set_yticks([200 * i for i in range(11)])
    # axes[1].set_yticklabels([2 * i for i in range(11)], fontsize=12)
    # axes[1].set_ylim(10, 950)
    axes[1].set_ylabel(r"$N_{contacts} (\times 100)$", fontsize=16)

    figname = "TDP43_contact_count_{0}.svg".format(phase_state)
    plt.savefig(figname, dpi=150)

if __name__ == '__main__':
    main("dense")
    main("dilute")
