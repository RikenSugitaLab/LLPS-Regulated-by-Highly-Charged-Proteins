#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():
    n_T = 1
    n_C = 9
    n_run = 5

    fT = lambda i : i * 10 + 290
    fC = lambda i : i * 10 + 10

    fig, axes = plt.subplots(1, 1, figsize=(5, 5), constrained_layout=True, sharex=False, sharey=False)
    color_list = ["r", "g", "b", "m"]

    for i_r in range(n_run):
        i_T = 0
        c_tmp = color_list[i_T]
        lo_d1_list = []
        hi_d1_list = []
        lo_d2_list = []
        hi_d2_list = []
        for i_C in range(n_C):
            # if i_C % 2 == 0:
                # continue
            print("analyzing #T: ", i_T + 1, "  #Hero: ", i_C + 1)
            data_1_fname = "../01_density_in_box/tdp43_100_hero11_{0}_T{1}_r{2:>02d}_tdp43_density_averaged.dat".format(fC(i_C), fT(i_T), i_r+1)
            data_2_fname = "../01_density_in_box/tdp43_100_hero11_{0}_T{1}_r{2:>02d}_hero11_density_averaged.dat".format(fC(i_C), fT(i_T), i_r+1)

            # calculate average of density of TDP43
            data_1 = np.loadtxt(data_1_fname, usecols=(3, 4))[6000:, :]
            lo_d1_ave, hi_d1_ave = np.mean(data_1, axis=0)
            lo_d1_std, hi_d1_std = np.std(data_1, axis=0)


            # calculate average of density of Hero11
            data_2 = np.loadtxt(data_2_fname, usecols=(3, 4))[6000:, :]
            lo_d2_ave, hi_d2_ave = np.mean(data_2, axis=0)
            lo_d2_std, hi_d2_std = np.std(data_2, axis=0)

            axes.errorbar([lo_d1_ave], [lo_d2_ave], xerr=[lo_d1_std], yerr=[lo_d2_std], elinewidth=1.0, ecolor=c_tmp, capsize=3, capthick=1.0, marker="o", mfc="white", mec=c_tmp, mew=1.5, ms=5)
            axes.errorbar([hi_d1_ave], [hi_d2_ave], xerr=[hi_d1_std], yerr=[hi_d2_std], elinewidth=1.0, ecolor=c_tmp, capsize=3, capthick=1.0, marker="o", mfc="white", mec=c_tmp, mew=1.5, ms=5)
            axes.plot([lo_d1_ave, hi_d1_ave], [lo_d2_ave, hi_d2_ave], ls="--", c=c_tmp, alpha=0.3, lw=0.5)

            print(" >>>> ", hi_d1_ave * 154 + hi_d2_ave * 99)

            lo_d1_list.append(lo_d1_ave)
            hi_d1_list.append(hi_d1_ave)
            lo_d2_list.append(lo_d2_ave)
            hi_d2_list.append(hi_d2_ave)

        # ================
        # add 100:100 data
        # ================
        for i_C in [10]:
            print("analyzing #T: ", i_T + 1, "  #Hero: ", i_C)
            data_1_fname = "../../../../20220210_first_test_with_Hero11_helix/34_tdp43_100_hero11_wt_100_multi_T_250_350_hokusai_md4/analysis/density_evolution/tdp43_md4_T_{0:>02d}_density_averaged.dat".format(i_T + 2)
            data_2_fname = "../../../../20220210_first_test_with_Hero11_helix/34_tdp43_100_hero11_wt_100_multi_T_250_350_hokusai_md4/analysis/density_evolution/hero11_md4_T_{0:>02d}_density_averaged.dat".format(i_T + 2)

            # calculate average of density of TDP43
            data_1 = np.loadtxt(data_1_fname, usecols=(3, 4))
            lo_d1_ave, hi_d1_ave = np.mean(data_1, axis=0)
            lo_d1_std, hi_d1_std = np.std(data_1, axis=0)


            # calculate average of density of Hero11
            data_2 = np.loadtxt(data_2_fname, usecols=(3, 4))
            lo_d2_ave, hi_d2_ave = np.mean(data_2, axis=0)
            lo_d2_std, hi_d2_std = np.std(data_2, axis=0)

            axes.errorbar([lo_d1_ave], [lo_d2_ave], xerr=[lo_d1_std], yerr=[lo_d2_std], elinewidth=1.0, ecolor=c_tmp, capsize=3, capthick=1.0, marker="o", mfc="white", mec=c_tmp, mew=1.5, ms=5)
            axes.errorbar([hi_d1_ave], [hi_d2_ave], xerr=[hi_d1_std], yerr=[hi_d2_std], elinewidth=1.0, ecolor=c_tmp, capsize=3, capthick=1.0, marker="o", mfc="white", mec=c_tmp, mew=1.5, ms=5)
            axes.plot([lo_d1_ave, hi_d1_ave], [lo_d2_ave, hi_d2_ave], ls="--", c=c_tmp, alpha=0.3, lw=0.5)

            lo_d1_list.append(lo_d1_ave)
            hi_d1_list.append(hi_d1_ave)
            lo_d2_list.append(lo_d2_ave)
            hi_d2_list.append(hi_d2_ave)


        axes.plot(lo_d1_list, lo_d2_list, ls="-", c=c_tmp, lw=1.2, alpha = 0.5)
        axes.plot(hi_d1_list, hi_d2_list, ls="-", c=c_tmp, lw=1.2, alpha = 0.5)

        axes.set_xticks([0.005 * i for i in range(4)])
        axes.set_xticklabels([0, 5, 10, 15], fontsize=12)
        axes.set_xlim(-0.0005, 0.0155)
        axes.set_xlabel("density of TDP43 (mM)", fontsize=16)

        axes.set_yticks([0.001 * i for i in range(6)])
        axes.set_yticklabels([0, 1, 2, 3, 4, 5], fontsize=12)
        axes.set_ylim(-0.0002, 0.0048)
        axes.set_ylabel("density of Hero11 (mM)", fontsize=16)


    # =======
    # test...
    # =======
    x_test = np.linspace(0, 0.015)
    y1_test = -1.5 * (x_test - 0.0138)
    y2_test = -1.5 * (x_test - 0.0129)
    y3_test = -1.5 * (x_test - 0.0116)
    y4_test = -1.5 * (x_test - 0.0098)
    # axes.plot(x_test, y1_test, c="gray", ls="-.", lw=1.0, alpha=0.8)
    # axes.plot(x_test, y2_test, c="gray", ls="-.", lw=1.0, alpha=0.8)
    # axes.plot(x_test, y3_test, c="gray", ls="-.", lw=1.0, alpha=0.8)
    # axes.plot(x_test, y4_test, c="gray", ls="-.", lw=1.0, alpha=0.8)

    # add homotypic TDP43 data...
    homo_tdp43_data = np.loadtxt("../../../../../TDP43/20220206_first_test/25_tdp43_100_multi_temperature_260_350_HPS-Urry_FIX_BONDLENGTH/analysis/phase-diagram_new/phase_diagram_raw.dat", usecols=(1, 2))
    homo_tdp43_lo_density = homo_tdp43_data[:, 0]
    homo_tdp43_hi_density = homo_tdp43_data[:, 1]
    for i_T in range(n_T):
        axes.scatter(homo_tdp43_lo_density[i_T], [0], marker="P", c=color_list[i_T], s=50)
        axes.scatter(homo_tdp43_hi_density[i_T], [0], marker="x", c=color_list[i_T], s=50)

    # plt.show()
    # plt.savefig("TDP43_Hero11_binary_phase_diagram.svg", dpi=150)


if __name__ == '__main__':
    main()
