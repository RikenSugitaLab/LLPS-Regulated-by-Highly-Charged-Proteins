#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def plot_complex_contacts(n_hero, T, phase_state, i_run):

    name_pro1 = "tdp43"
    name_pro2 = "hero11"
    system_name = "{0}_{1}_{2}_{3}".format(name_pro1, 100, name_pro2, n_hero)

    print("Plotting contact map for system: ", system_name, " at T = ", T, " in phase: ", phase_state)

    n_res_pro1 = 154
    n_res_pro2 = 99

    # ========================
    # load 2d contact matrices
    # ========================
    # p1-p2
    p1_p2_2d_fname = "./md{0}/{1}_t{2:3d}_md{0}.{3}-{4}.{5}.contact_matrix.dat".format(i_run, system_name, T, name_pro1, name_pro2, phase_state)
    cm_12 = np.loadtxt(p1_p2_2d_fname)
    cm_12_zmax = np.max(cm_12)
    # cm_12_vmax = round(cm_12_zmax, int(-np.log10(cm_12_zmax)) + 1) if cm_12_zmax > 1e-6 else 1e-5
    cm_12_vmax = 0.05
    # p1-p1
    p1_p1_2d_fname = "./md{0}/{1}_t{2:3d}_md{0}.{3}-{4}.{5}.contact_matrix.dat".format(i_run, system_name, T, name_pro1, name_pro1, phase_state)
    cm_11 = np.loadtxt(p1_p1_2d_fname)
    cm_11_zmax = np.max(cm_11)
    # cm_11_vmax = round(cm_11_zmax, int(-np.log10(cm_11_zmax)) + 1) if cm_11_zmax > 1e-6 else 1e-5
    cm_11_vmax = 0.25
    # p2-p2
    p2_p2_2d_fname = "./md{0}/{1}_t{2:3d}_md{0}.{3}-{4}.{5}.contact_matrix.dat".format(i_run, system_name, T, name_pro2, name_pro2, phase_state)
    cm_22 = np.loadtxt(p2_p2_2d_fname)
    cm_22_zmax = np.max(cm_22)
    # cm_22_vmax = round(cm_22_zmax, int(-np.log10(cm_22_zmax)) + 1) if cm_22_zmax > 1e-6 else 1e-5
    cm_22_vmax = 0.05

    # ======================
    # load 1d contact counts
    # ======================
    # p1-p2
    p1_p2_1d_fname = "./md{0}_1d_count/{1}_1d_contact_count.{1}-{2}.{3}.T{4}.{5}_{6}.dat".format(i_run, name_pro1, name_pro2, phase_state, T, name_pro2, n_hero)
    cc_12 = np.loadtxt(p1_p2_1d_fname)
    # p2-p1
    p2_p1_1d_fname = "./md{0}_1d_count/{1}_1d_contact_count.{1}-{2}.{3}.T{4}.{5}_{6}.dat".format(i_run, name_pro2, name_pro1, phase_state, T, name_pro2, n_hero)
    cc_21 = np.loadtxt(p2_p1_1d_fname)
    # p1-p1
    p1_p1_1d_fname = "./md{0}_1d_count/{1}_1d_contact_count.{1}-{2}.{3}.T{4}.{5}_{6}.dat".format(i_run, name_pro1, name_pro1, phase_state, T, name_pro2, n_hero)
    cc_11 = np.loadtxt(p1_p1_1d_fname)
    # p2-p2
    p2_p2_1d_fname = "./md{0}_1d_count/{1}_1d_contact_count.{1}-{2}.{3}.T{4}.{5}_{6}.dat".format(i_run, name_pro2, name_pro2, phase_state, T, name_pro2, n_hero)
    cc_22 = np.loadtxt(p2_p2_1d_fname)


    ###########################################################################
    #                                 plotting                                #
    ###########################################################################
    fig = plt.figure(figsize=(9, 9))

    # ==================
    # prepare rectangles
    # ==================
    cm_pixel = 0.003
    width = n_res_pro1 * cm_pixel
    height = n_res_pro2 * cm_pixel
    left = 0.1
    bottom = 0.1
    spacing = 0.01

    rect_cm12 = [left, bottom, width, height]
    rect_cm11 = [left, bottom + height + spacing, width, width]
    rect_cm22 = [left + width + spacing, bottom, height, height]
    rect_cc1  = [left, bottom + height + width + 2 * spacing, width, 0.1]
    rect_cc2  = [left + height + width + 2 * spacing, bottom, 0.1, height]

    # =======================
    # contact map pro1 - pro2
    # =======================
    ax_cm_12 = plt.axes(rect_cm12, facecolor='w')
    # x-axis
    ax_cm_12.set_xticks([20 * i - 1 for i in range(8)])
    ax_cm_12.set_xticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    ax_cm_12.set_xlim(-0.5, 153.5)
    ax_cm_12.set_xlabel("residue index (TDP43)", fontsize=16)
    # y-axis
    ax_cm_12.set_yticks([20 * i - 1 for i in range(8)])
    ax_cm_12.set_yticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    ax_cm_12.set_ylim(-0.5, 98.5)
    ax_cm_12.set_ylabel("residue index (Hero11)", fontsize=16)
    # plot contact map
    ax_cm_12.imshow(np.transpose(cm_12), cmap="Purples", vmin=1e-5, vmax=cm_12_vmax, origin="lower", aspect="equal")

    # =======================
    # contact map pro1 - pro1
    # =======================
    ax_cm_11 = plt.axes(rect_cm11, facecolor='w', sharex=ax_cm_12)
    # x-axis
    # ax_cm_11.set_xticks([20 * i - 1 for i in range(8)])
    # ax_cm_11.set_xticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    # ax_cm_11.set_xlim(-0.5, 153.5)
    # ax_cm_11.set_xlabel("residue index (TDP43)", fontsize=16)
    # y-axis
    ax_cm_11.set_yticks([20 * i - 1 for i in range(8)])
    ax_cm_11.set_yticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    ax_cm_11.set_ylim(-0.5, 153.5)
    ax_cm_11.set_ylabel("residue index (TDP43)", fontsize=16)
    # plot contact map
    ax_cm_11.imshow(cm_11, cmap="Blues", vmin=1e-5, vmax=cm_11_vmax, origin="lower", aspect="equal")
    # remove bottom label
    ax_cm_11.tick_params(axis="x", labelbottom=False)

    # =======================
    # contact map pro2 - pro2
    # =======================
    ax_cm_22 = plt.axes(rect_cm22, facecolor='w', sharey=ax_cm_12)
    # x-axis
    ax_cm_22.set_xticks([20 * i - 1 for i in range(8)])
    ax_cm_22.set_xticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    ax_cm_22.set_xlim(-0.5, 98.5)
    ax_cm_22.set_xlabel("residue index (Hero11)", fontsize=16)
    # y-axis
    # ax_cm_22.set_yticks([20 * i - 1 for i in range(8)])
    # ax_cm_22.set_yticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    # ax_cm_22.set_ylim(-0.5, 98.5)
    # ax_cm_22.set_ylabel("residue index (Hero11)", fontsize=16)
    # plot contact map
    ax_cm_22.imshow(cm_22, cmap="Reds", vmin=1e-5, vmax=cm_22_vmax, origin="lower", aspect="equal")
    # remove left label
    ax_cm_22.tick_params(axis="y", labelleft=False)


    # =====================
    # contact count of pro1
    # =====================
    # add the 1d contact number count axis
    ax_cc_1 = plt.axes(rect_cc1, facecolor='w', sharex=ax_cm_12)
    # plot 1d contact probabilities
    X = cc_11[:, 0] - 1
    Y11 = cc_11[:, 1]
    ax_cc_1.plot(X, Y11, color="b", ls="--", lw=1.0)
    Y12 = cc_12[:, 1]
    ax_cc_1.plot(X, Y12, color="r", ls="--", lw=1.0)
    Y_all = cc_11[:, 1] + cc_12[:, 1]
    ax_cc_1.plot(X, Y_all, color="b", ls="-", lw=2.0)
    # x-axis
    # ax_cc_1.set_xticks([20 * i - 1 for i in range(8)])
    # ax_cc_1.set_xticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    # ax_cc_1.set_xlim(-0.5, 153.5)
    # ax_cc_1.set_xlabel("residue index (TDP43)", fontsize=16)
    # y-axis
    ax_cc_1.set_yticks([20 * i - 1 for i in range(8)])
    ax_cc_1.set_yticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    ax_cc_1.set_ylim(-0.5, 40.5)
    ax_cc_1.set_ylabel("N (TDP43)", fontsize=16)
    # remove bottom label
    ax_cc_1.tick_params(axis="x", labelbottom=False)

    # =====================
    # contact count of pro1
    # =====================
    # add the 1d contact number count axis
    ax_cc_2 = plt.axes(rect_cc2, facecolor='w', sharey=ax_cm_12)
    # plot 1d contact probabilities
    Y = cc_22[:, 0] - 1
    X22 = cc_22[:, 1]
    ax_cc_2.plot(X22, Y, color="r", ls="--", lw=1.0)
    X21 = cc_21[:, 1]
    ax_cc_2.plot(X21, Y, color="b", ls="--", lw=1.0)
    X_all = cc_22[:, 1] + cc_21[:, 1]
    ax_cc_2.plot(X_all, Y, color="r", ls="-", lw=2.0)
    # x-axis
    ax_cc_2.set_xticks([20 * i - 1 for i in range(8)])
    ax_cc_2.set_xticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    ax_cc_2.set_xlim(-0.5, 40.5)
    ax_cc_2.set_xlabel("N (Hero11)", fontsize=16)
    # y-axis
    # ax_cc_2.set_yticks([20 * i - 1 for i in range(8)])
    # ax_cc_2.set_yticklabels(["", 20, 40, 60, 80, 100, 120, 140], fontsize=10)
    # ax_cc_2.set_ylim(-0.5, 153.5)
    # ax_cc_2.set_ylabel("residue index (Hero11)", fontsize=16)
    # remove bottom label
    ax_cc_2.tick_params(axis="y", labelleft=False)

    # output some important info
    ax_cm_12.text(160, 120, "cm_12_vmax = {0:8.3f}".format(cm_12_vmax))
    ax_cm_12.text(160, 140, "cm_11_vmax = {0:8.3f}".format(cm_11_vmax))
    ax_cm_12.text(160, 160, "cm_22_vmax = {0:8.3f}".format(cm_22_vmax))

    # plt.show()
    figname = "md{0}_T{1}_{2}_{3}_{4}.contact_plot.svg".format(i_run, T,  phase_state, name_pro2, n_hero)
    plt.savefig(figname, dpi=150)

if __name__ == '__main__':
    for i_run in [4]:
    # for i_run in range(1):
        for i_t in [0]:
        # for i_t in range(1):
            plot_complex_contacts(i_run * 10 + 10, i_t * 10 + 260, "dense", 3)
            plot_complex_contacts(i_run * 10 + 10, i_t * 10 + 260, "dilute", 3)
