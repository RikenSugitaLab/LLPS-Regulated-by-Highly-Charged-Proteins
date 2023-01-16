#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():

    T_list = [290, 295]
    n_list = [10 * i + 10 for i in range(9)]

    pro1_name = "tdp43"
    pro2_name = "hero11"
    n_pro1 = 100

    fig, axes = plt.subplots(2, 1, figsize=(5, 5), constrained_layout=True, sharex=True, sharey=False)
    color_list = ["r", "g", "b", "m"]
    rg_0_list = [31.126, 31.126]
    for i_T, T in enumerate(T_list):

        of1_dense  = open("{0}_T{1}.rg.dense.dat".format( pro1_name, T), "w")
        of2_dense  = open("{0}_T{1}.rg.dense.dat".format( pro2_name, T), "w")
        of1_dilute = open("{0}_T{1}.rg.dilute.dat".format(pro1_name, T), "w")
        of2_dilute = open("{0}_T{1}.rg.dilute.dat".format(pro2_name, T), "w")

        y_p1_rg_dense_ave = []
        y_p2_rg_dense_ave = []
        y_p1_rg_dense_err = []
        y_p2_rg_dense_err = []
        y_p1_rg_dilute_ave = []
        y_p2_rg_dilute_ave = []
        y_p1_rg_dilute_err = []
        y_p2_rg_dilute_err = []

        for n in n_list:
            system_name = "{0}_{1}_{2}_{3}".format(pro1_name, n_pro1, pro2_name, n)

            print(" Dealing with system: ", system_name)

            pro1_rg_dense_all = np.array([])
            pro2_rg_dense_all = np.array([])
            pro1_rg_dilute_all = np.array([])
            pro2_rg_dilute_all = np.array([])

            for irun in [1,2,3,4,5]:
                for j in [2, 3]:
                    # load rg data
                    rg_data_fname = "./md{0}/{1}_T{2}_md{0}_r{3:0>2d}.{4}_rg.dat"
                    print(" >> Processing ", irun, " - ", j)
                    pro1_rg = np.loadtxt(rg_data_fname.format(j, system_name, T_list[i_T], irun, pro1_name))
                    pro2_rg = np.loadtxt(rg_data_fname.format(j, system_name, T_list[i_T], irun, pro2_name))

                    # load position data
                    pos_data_fname = "../03_molecular_position/md{0}/{1}_T{2}_md{0}_r{3:0>2d}.{4}_pos.dat"
                    pro1_pos = np.loadtxt(pos_data_fname.format(j, system_name, T_list[i_T], irun, pro1_name))
                    pro2_pos = np.loadtxt(pos_data_fname.format(j, system_name, T_list[i_T], irun, pro2_name))

                    # pick out molecules in condensate
                    pro1_rg_dense = pro1_rg[pro1_pos > 0.5]
                    pro2_rg_dense = pro2_rg[pro2_pos > 0.5]
                    pro1_rg_dense_all = np.concatenate((pro1_rg_dense_all, pro1_rg_dense))
                    pro2_rg_dense_all = np.concatenate((pro2_rg_dense_all, pro2_rg_dense))

                    # pick out molecules outside condensate
                    pro1_rg_dilute = pro1_rg[pro1_pos < 0.5]
                    pro2_rg_dilute = pro2_rg[pro2_pos < 0.5]
                    pro1_rg_dilute_all = np.concatenate((pro1_rg_dilute_all, pro1_rg_dilute))
                    pro2_rg_dilute_all = np.concatenate((pro2_rg_dilute_all, pro2_rg_dilute))

            # check data shape...
            print(pro1_rg_dense_all.shape)
            print(pro2_rg_dense_all.shape)
            print(pro1_rg_dilute_all.shape)
            print(pro2_rg_dilute_all.shape)

            # ave and std
            p1_rg_dense_ave  = np.mean(pro1_rg_dense_all)
            p2_rg_dense_ave  = np.mean(pro2_rg_dense_all)
            p1_rg_dilute_ave = np.mean(pro1_rg_dilute_all)
            p2_rg_dilute_ave = np.mean(pro2_rg_dilute_all)

            p1_rg_dense_std  = np.std(pro1_rg_dense_all)
            p2_rg_dense_std  = np.std(pro2_rg_dense_all)
            p1_rg_dilute_std = np.std(pro1_rg_dilute_all)
            p2_rg_dilute_std = np.std(pro2_rg_dilute_all)

            p1_rg_dense_err  = p1_rg_dense_std  / np.size(pro1_rg_dense_all)**0.5
            p2_rg_dense_err  = p2_rg_dense_std  / np.size(pro2_rg_dense_all)**0.5
            p1_rg_dilute_err = p1_rg_dilute_std / np.size(pro1_rg_dilute_all)**0.5
            p2_rg_dilute_err = p2_rg_dilute_std / np.size(pro2_rg_dilute_all)**0.5

            # output data
            of1_dense.write( "{0:3d}  {1:8.3f}   {2:8.3f}   {3:8.3f} \n".format(n, p1_rg_dense_ave, p1_rg_dense_std, p1_rg_dense_err ))
            of2_dense.write( "{0:3d}  {1:8.3f}   {2:8.3f}   {3:8.3f} \n".format(n, p2_rg_dense_ave, p2_rg_dense_std, p2_rg_dense_err ))
            of1_dilute.write("{0:3d}  {1:8.3f}   {2:8.3f}   {3:8.3f} \n".format(n, p1_rg_dilute_ave, p1_rg_dilute_std, p1_rg_dilute_err ))
            of2_dilute.write("{0:3d}  {1:8.3f}   {2:8.3f}   {3:8.3f} \n".format(n, p2_rg_dilute_ave, p2_rg_dilute_std, p2_rg_dilute_err ))

            y_p1_rg_dense_ave.append(p1_rg_dense_ave)
            y_p2_rg_dense_ave.append(p2_rg_dense_ave)
            y_p1_rg_dense_err.append(p1_rg_dense_err)
            y_p2_rg_dense_err.append(p2_rg_dense_err)
            y_p1_rg_dilute_ave.append(p1_rg_dilute_ave)
            y_p2_rg_dilute_ave.append(p2_rg_dilute_ave)
            y_p1_rg_dilute_err.append(p1_rg_dilute_err)
            y_p2_rg_dilute_err.append(p2_rg_dilute_err)

        of1_dense.close()
        of2_dense.close()
        of1_dilute.close()
        of2_dilute.close()


        # PLOT!!!
        axes[0].errorbar(n_list, y_p1_rg_dense_ave, yerr=y_p1_rg_dense_err, ls="-", c=color_list[i_T], elinewidth=1.5, ecolor=color_list[i_T], capsize=5, capthick=1.5, marker="o", mfc="none", mec=color_list[i_T], mew=1.5, ms=5)
        axes[0].errorbar(n_list, y_p1_rg_dilute_ave, yerr=y_p1_rg_dilute_err, ls="--", c=color_list[i_T], elinewidth=1.5, ecolor=color_list[i_T], capsize=5, capthick=1.5, marker="o", mfc="none", mec=color_list[i_T], mew=1.5, ms=5)

        axes[0].axhline(y=rg_0_list[i_T], ls=":", c=color_list[i_T])

        axes[1].errorbar(n_list, y_p2_rg_dense_ave, yerr=y_p2_rg_dense_err, ls="-", c=color_list[i_T], elinewidth=1.5, ecolor=color_list[i_T], capsize=5, capthick=1.5, marker="s", mfc="none", mec=color_list[i_T], mew=1.5, ms=5)
        axes[1].errorbar(n_list, y_p2_rg_dilute_ave, yerr=y_p2_rg_dilute_err, ls="--", c=color_list[i_T], elinewidth=1.5, ecolor=color_list[i_T], capsize=5, capthick=1.5, marker="s", mfc="none", mec=color_list[i_T], mew=1.5, ms=5)

    # axis styles
    axes[0].set_xticks(n_list)
    # axes[0].set_xticklabels(n_list, fontsize=12)
    axes[0].set_xlim(5, 95)
    # axes[0].set_xlabel(r"$n_{Hero11}$", fontsize=16)
    axes[0].set_yticks([28, 30, 32, 34])
    axes[0].set_yticklabels(["28.0", "30.0", "32.0", "34.0"], fontsize=12)
    axes[0].set_ylim(27.9, 34.5)
    axes[0].set_ylabel(r"$R_g(TDP43) (\AA)$", fontsize=16)

    axes[1].set_xticks(n_list)
    axes[1].set_xticklabels(n_list, fontsize=12)
    axes[1].set_xlim(5, 95)
    axes[1].set_xlabel(r"$n_{Hero11}$", fontsize=16)
    axes[1].set_yticks([33, 33.5, 34])
    axes[1].set_yticklabels(["33.0", "33.5", "34.0"], fontsize=12)
    axes[1].set_ylim(32.6, 34.4)
    axes[1].set_ylabel(r"$R_g(TDP43) (\AA)$", fontsize=16)

    fig_name = "{0}_Rg_all_vs_n_hero.svg".format(system_name)
    plt.savefig(fig_name, dpi=150)

if __name__ == '__main__':
    main()
