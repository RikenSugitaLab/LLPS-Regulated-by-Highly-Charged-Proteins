#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main(n_hero):

    print("Analyzing systems with #Hero11 = ", n_hero)

    # params
    n_temperature = 4
    T_list = [i * 10 + 260 for i in range(n_temperature)]
    n_run = 3

    n_bins = 25

    pro1_name = "tdp43"
    pro2_name = "hero11"
    n_pro1 = 100
    n_pro2 = n_hero

    system_name = "{0}_{1}_{2}_{3}".format(pro1_name, n_pro1, pro2_name, n_pro2)


    # prepare the plotting ...
    fig, axes = plt.subplots(2, n_temperature, figsize=(4 * n_temperature, 5), constrained_layout=True, sharex=True, sharey=True)

    # =======================================================
    # calculate rg distributions based on molecular positions
    # =======================================================
    for i in range(n_temperature):
        pro1_rg_dense_all = np.array([])
        pro2_rg_dense_all = np.array([])
        pro1_rg_dilute_all = np.array([])
        pro2_rg_dilute_all = np.array([])
        for j in [2, 3]:
            # load rg data
            rg_data_fname = "./md{0}/{1}_t{2}_md{0}.{3}_rg.dat"
            pro1_rg = np.loadtxt(rg_data_fname.format(j, system_name, T_list[i], pro1_name))
            pro2_rg = np.loadtxt(rg_data_fname.format(j, system_name, T_list[i], pro2_name))

            # load position data
            pos_data_fname = "../03_molecular_position/md{0}/{1}_t{2}_md{0}.{3}_pos.dat"
            pro1_pos = np.loadtxt(pos_data_fname.format(j, system_name, T_list[i], pro1_name))
            pro2_pos = np.loadtxt(pos_data_fname.format(j, system_name, T_list[i], pro2_name))

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
        # print(pro1_rg_dense_all.shape)
        # print(pro2_rg_dense_all.shape)
        # print(pro1_rg_dilute_all.shape)
        # print(pro2_rg_dilute_all.shape)

        # histogram
        p1_rg_dense_hist, p1_rg_dense_edge = np.histogram(pro1_rg_dense_all, bins=n_bins, density=True)
        p2_rg_dense_hist, p2_rg_dense_edge = np.histogram(pro2_rg_dense_all, bins=n_bins, density=True)
        p1_rg_dilute_hist, p1_rg_dilute_edge = np.histogram(pro1_rg_dilute_all, bins=n_bins, density=True)
        p2_rg_dilute_hist, p2_rg_dilute_edge = np.histogram(pro2_rg_dilute_all, bins=n_bins, density=True)

        # output histogram
        of1_dense  = open("{0}_T{1}.{2}.rg_distribution_dense.dat".format(system_name, T_list[i], pro1_name), "w")
        of2_dense  = open("{0}_T{1}.{2}.rg_distribution_dense.dat".format(system_name, T_list[i], pro2_name), "w")
        of1_dilute = open("{0}_T{1}.{2}.rg_distribution_dilute.dat".format(system_name, T_list[i], pro1_name), "w")
        of2_dilute = open("{0}_T{1}.{2}.rg_distribution_dilute.dat".format(system_name, T_list[i], pro2_name), "w")
        for k in range(n_bins):
            of1_dense.write("{0:8.3f}   {1:8.3f} \n".format( 0.5*(p1_rg_dense_edge[k] + p1_rg_dense_edge[k + 1]), p1_rg_dense_hist[k]))
            of2_dense.write("{0:8.3f}   {1:8.3f} \n".format( 0.5*(p2_rg_dense_edge[k] + p2_rg_dense_edge[k + 1]), p2_rg_dense_hist[k]))
            of1_dilute.write("{0:8.3f}   {1:8.3f} \n".format( 0.5*(p1_rg_dilute_edge[k] + p1_rg_dilute_edge[k + 1]), p1_rg_dilute_hist[k]))
            of2_dilute.write("{0:8.3f}   {1:8.3f} \n".format( 0.5*(p2_rg_dilute_edge[k] + p2_rg_dilute_edge[k + 1]), p2_rg_dilute_hist[k]))
        of1_dense.close()
        of2_dense.close()
        of1_dilute.close()
        of2_dilute.close()


        # PLOT!!!
        X1_dense = 0.5 * (p1_rg_dense_edge[1:] + p1_rg_dense_edge[:-1])
        X2_dense = 0.5 * (p2_rg_dense_edge[1:] + p2_rg_dense_edge[:-1])
        X1_dilute = 0.5 * (p1_rg_dilute_edge[1:] + p1_rg_dilute_edge[:-1])
        X2_dilute = 0.5 * (p2_rg_dilute_edge[1:] + p2_rg_dilute_edge[:-1])

        axes[0, i].plot(X1_dense, p1_rg_dense_hist, "b-", lw=2)
        axes[0, i].plot(X1_dilute, p1_rg_dilute_hist, "b--", lw=1.5)
        axes[1, i].plot(X2_dense, p2_rg_dense_hist, "r-", lw=2)
        axes[1, i].plot(X2_dilute, p2_rg_dilute_hist, "r--", lw=1.5)

        axes[0, i].set_xticks([10 * k for k in range(11)])
        axes[0, i].set_xlim(12, 78)
        axes[1, i].set_xticks([10 * k for k in range(11)])
        axes[1, i].set_xlim(12, 78)
        axes[1, i].set_xticklabels([10 * k for k in range(11)], fontsize=12)
        axes[1, i].set_xlabel(r"$R_g (\AA)$", fontsize=16)

        axes[0, i].set_yticks([0.02 * k for k in range(5)])
        axes[0, i].set_ylim(0, 0.078)
        axes[1, i].set_yticks([0.02 * k for k in range(5)])
        axes[1, i].set_ylim(0, 0.078)
        if i == 0:
            axes[0, i].set_yticklabels([0, "0.02", "0.04", "0.06", "0.08"])
            axes[1, i].set_yticklabels([0, "0.02", "0.04", "0.06", "0.08"])
            axes[0, i].set_ylabel("probability", fontsize=16)
            axes[1, i].set_ylabel("probability", fontsize=16)

    fig_name = "{0}_Rg_distributions_all_temperatures.svg".format(system_name)
    plt.savefig(fig_name, dpi=150)

if __name__ == '__main__':
    for i in range(9):
        main(i * 10 + 10)
