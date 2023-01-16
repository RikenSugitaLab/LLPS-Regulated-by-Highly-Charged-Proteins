#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():
    num_temperature = 6
    num_temperature_meaningful = 3
    num_run = 5
    MD_T_list = [120 + i * 5 for i in range(num_temperature)]
    meaningful_T_list = [i for i in range(num_temperature_meaningful)]

    # ==================================
    # read in temperatures and densities
    # ==================================
    ana_T_list = np.zeros(num_temperature_meaningful)
    ana_lo_density = np.zeros(num_temperature_meaningful)
    ana_hi_density = np.zeros(num_temperature_meaningful)
    ana_delta_density = np.zeros(num_temperature_meaningful)

    for j in range(num_temperature_meaningful):
        for k in range(num_run):
            fname = "../density_evolution/hero11_alpha_100_md1_T_{0:0>2d}_{1:0>2d}_density_averaged.dat".format(j + 1, k + 1)
            if k == 0:
                density_data_all = np.loadtxt(fname, usecols=(3, 4))[500:, :]
            else:
                density_data = np.loadtxt(fname, usecols=(3, 4))[500:, :]
                density_data_all = np.append(density_data_all, density_data, axis=0)
        lo_d_ave, hi_d_ave = np.mean(density_data_all, axis=0)
        ana_T_list[j] = MD_T_list[j]
        ana_lo_density[j] = lo_d_ave
        ana_hi_density[j] = hi_d_ave
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
    ax.scatter(ana_lo_density, ana_T_list, marker="o", s=40, c="black")
    ax.scatter(ana_hi_density, ana_T_list, marker="o", s=40, c="black")
    ax.scatter([RHO_C], [T_C], marker="^", s=40, c="gray")

    ax.set_xticks([0.005 * k for k in range(10)])
    ax.set_xticklabels([0.005 * k for k in range(10)], fontsize=12)
    ax.set_xlim(-0.001, 0.047)
    ax.set_xlabel("density (M)", fontsize=16)
    ax.set_yticks(MD_T_list)
    ax.set_yticklabels(MD_T_list, fontsize=12)
    ax.set_ylim(T_min, T_C * 1.02)
    ax.set_ylabel("temperature (K)", fontsize=16)

    figname = "phase_diagram_fitted.svg"
    plt.savefig(figname)
    # plt.show()

    # ===========================
    # output data for future plot
    # ===========================
    raw_data_merge = np.vstack((ana_T_list, ana_lo_density, ana_hi_density))
    fit_data_merge = np.vstack((fit_T_list, fit_lo_d_list,  fit_hi_d_list))

    np.savetxt("phase_diagram_raw.dat", raw_data_merge.T)
    np.savetxt("phase_diagram_fit.dat", fit_data_merge.T)


if __name__ == '__main__':
    main()
