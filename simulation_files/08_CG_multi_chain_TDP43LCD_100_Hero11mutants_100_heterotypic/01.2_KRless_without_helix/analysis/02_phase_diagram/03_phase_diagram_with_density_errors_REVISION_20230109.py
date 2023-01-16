#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():
    pro0_name = "tdp43"
    num_temperature = 10
    num_temperature_meaningful = 4
    MD_T_list = [260 + i * 10 for i in range(num_temperature)]
    # MD_T_list = [250 + i * 10 for i in range(num_temperature)]
    meaningful_T_list = [i for i in range(num_temperature_meaningful)]

    # ==================================
    # read in temperatures and densities
    # ==================================
    ana_T_list = np.zeros(num_temperature_meaningful)
    ana_lo_density = np.zeros(num_temperature_meaningful)
    ana_hi_density = np.zeros(num_temperature_meaningful)
    ana_lo_std = np.zeros(num_temperature_meaningful)
    ana_hi_std = np.zeros(num_temperature_meaningful)
    ana_lo_err = np.zeros(num_temperature_meaningful)
    ana_hi_err = np.zeros(num_temperature_meaningful)
    ana_delta_density = np.zeros(num_temperature_meaningful)

    for j in range(num_temperature_meaningful):
        fname = "../density_evolution/{0}_md2_T_{1:0>2d}_density_averaged.dat".format(pro0_name, j + 2)
        density_data = np.loadtxt(fname, usecols=(3, 4))
        lo_d_ave, hi_d_ave = np.mean(density_data, axis=0)
        lo_d_std, hi_d_std = np.std(density_data, axis=0)
        ana_T_list[j] = MD_T_list[j]
        ana_lo_density[j] = lo_d_ave
        ana_hi_density[j] = hi_d_ave
        ana_lo_std[j] = lo_d_std
        ana_hi_std[j] = hi_d_std
        print("data size: ", np.shape(density_data)[0])
        ana_lo_err[j] = lo_d_std / np.shape(density_data)[0]**0.5
        ana_hi_err[j] = hi_d_std / np.shape(density_data)[0]**0.5
        ana_delta_density[j] = hi_d_ave - lo_d_ave

    # fit 1: critical temperature
    X = ana_delta_density**(1/0.325)
    Y = ana_T_list
    tc_fit = np.polyfit(X, Y, 1)
    T_C = tc_fit[1]
    D_rho = (- T_C / tc_fit[0]) ** 0.325
    print(D_rho, T_C)

    # fit 2: critical density
    X = -1.0 * ana_T_list + T_C
    Y = 0.5 * (ana_lo_density + ana_hi_density)
    dc_fit = np.polyfit(X, Y, 1)
    A = dc_fit[0]
    RHO_C = dc_fit[1]
    print(A, RHO_C)

    # ==================
    # plot phase diagram
    # ==================
    fig, ax = plt.subplots(1, 1, figsize=(8, 5), constrained_layout=False, sharex=False, sharey=False)

    fit_num_points = 20
    T_min = MD_T_list[0] * 0.98
    g_space = np.geomspace(T_min / T_C, 1.0, num=100)
    g_space = (T_min / T_C) / g_space
    fit_T_list = T_C * g_space[::-1]
    fit_delta_dens = D_rho * ( 1.0 - fit_T_list / T_C)**0.325
    fit_sum_dens = 2 * (RHO_C + A * (T_C - fit_T_list))
    fit_lo_d_list = 0.5 * (fit_sum_dens - fit_delta_dens)
    fit_hi_d_list = 0.5 * (fit_delta_dens + fit_sum_dens)

    # plotting fitted curve
    ax.plot(fit_lo_d_list, fit_T_list, ls="-", lw=2.0, c="red", alpha=0.3)
    ax.plot(fit_hi_d_list, fit_T_list, ls="-", lw=2.0, c="red", alpha=0.3)

    # plotting raw data
    # ax.scatter(ana_lo_density, ana_T_list, marker="o", s=40, c="black")
    # ax.scatter(ana_hi_density, ana_T_list, marker="o", s=40, c="black")

    ax.errorbar(ana_lo_density, ana_T_list, xerr=ana_lo_std, elinewidth=1.0, capsize=3, capthick=1.0, marker="o", mfc="red", mec="red", mew=1.5, ms=5)
    # ax.errorbar(ana_lo_density, ana_T_list, xerr=ana_lo_err, elinewidth=1.0, capsize=5, capthick=1.0, marker="o", mfc="red", mec="red", mew=1.5, ms=5)
    ax.errorbar(ana_hi_density, ana_T_list, xerr=ana_hi_std, elinewidth=1.0, capsize=3, capthick=1.0, marker="o", mfc="red", mec="red", mew=1.5, ms=5)
    # ax.errorbar(ana_hi_density, ana_T_list, xerr=ana_hi_err, elinewidth=1.0, capsize=5, capthick=1.0, marker="o", mfc="red", mec="red", mew=1.5, ms=5)
    ax.scatter([RHO_C], [T_C], marker="^", s=40, c="gray")

    ax.set_xticks([0.005 * k for k in range(10)])
    ax.set_xticklabels([0.005 * k for k in range(10)], fontsize=12)
    ax.set_xlim(-0.001, 0.017)
    ax.set_xlabel("density (M)", fontsize=16)
    ax.set_yticks(MD_T_list)
    ax.set_yticklabels(MD_T_list, fontsize=12)
    # ax.set_ylim(T_min, T_C * 1.05)
    ax.set_ylim(255, 335)
    ax.set_ylabel("temperature (K)", fontsize=16)

    # figname = "{0}_phase_diagram_fitted.svg".format(pro0_name)
    figname = "20230109_REVISION_phase_diagram_fitted.svg"
    plt.savefig(figname)
    # plt.show()

    # ===========================
    # output data for future plot
    # ===========================
    # raw_data_merge = np.vstack((ana_T_list, ana_lo_density, ana_hi_density))
    # fit_data_merge = np.vstack((fit_T_list, fit_lo_d_list,  fit_hi_d_list))

    raw_data_merge = np.vstack((ana_T_list, ana_lo_density, ana_lo_std, ana_lo_err,  ana_hi_density, ana_hi_std, ana_hi_err))
    fit_data_merge = np.vstack((fit_T_list, fit_lo_d_list,  fit_hi_d_list))

    np.savetxt("20230109_REVISION_phase_diagram_raw.dat", raw_data_merge.T)
    np.savetxt("20230109_REVISION_phase_diagram_fit.dat", fit_data_merge.T)

    # pd_raw_fname = "{0}_phase_diagram_raw.dat".format(pro0_name)
    # pd_fit_fname = "{0}_phase_diagram_fit.dat".format(pro0_name)
    # np.savetxt(pd_raw_fname, raw_data_merge.T)
    # np.savetxt(pd_fit_fname, fit_data_merge.T)


if __name__ == '__main__':
    main()
