#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():
    pro0_name = "tdp43"
    pro1_name = "hero11"

    T_list = [250 + i * 10 for i in range(10)]
    num_temperatures = len(T_list)
    num_run = 1

    num_steps = 2500
    step_interval = 5
    plot_steps = num_steps // step_interval

    LOW_DENSITY_THRESHOLD = 0.003
    num_bins = 100


    fig, axes = plt.subplots(num_temperatures * num_run * 2, 1, figsize=(9, num_temperatures * num_run * 2), constrained_layout=True, sharex=True, sharey=False)
    for i, t in enumerate(T_list):
        print("Temperature:", t)
        fname0 = "{0}_md4_T_{1:0>2d}.density.dat".format(pro0_name, i + 1)
        fname1 = "{0}_md4_T_{1:0>2d}.density.dat".format(pro1_name, i + 1)
        data0_local = np.loadtxt(fname0)[::step_interval, :]
        data1_local = np.loadtxt(fname1)[::step_interval, :]

        pro0_density_profile_output_fname = "{0}_md4_T_{1:0>2d}_density_centered.dat".format(pro0_name, i + 1)
        pro0_density_average_output_fname = "{0}_md4_T_{1:0>2d}_density_averaged.dat".format(pro0_name, i + 1)
        pro1_density_profile_output_fname = "{0}_md4_T_{1:0>2d}_density_centered.dat".format(pro1_name, i + 1)
        pro1_density_average_output_fname = "{0}_md4_T_{1:0>2d}_density_averaged.dat".format(pro1_name, i + 1)
        pro0_dpf = open(pro0_density_profile_output_fname, "w")
        pro0_daf = open(pro0_density_average_output_fname, "w")
        pro1_dpf = open(pro1_density_profile_output_fname, "w")
        pro1_daf = open(pro1_density_average_output_fname, "w")

        # ================
        # guess the center
        # ================
        pro0_new_data = np.zeros((plot_steps, num_bins))
        pro1_new_data = np.zeros((plot_steps, num_bins))
        for j in range(plot_steps):
            # find the peak
            data0_j = data0_local[j, :]
            data1_j = data1_local[j, :]
            arr0_tmp = np.append(data0_j, data0_j)
            arr0_tmp = np.append(data0_j, arr0_tmp)
            arr1_tmp = np.append(data1_j, data1_j)
            arr1_tmp = np.append(data1_j, arr1_tmp)
            i_peak = np.argmax(data0_j) + num_bins
            i_lo_bound, i_hi_bound = i_peak, i_peak
            while arr0_tmp[i_lo_bound] >= LOW_DENSITY_THRESHOLD:
                i_lo_bound -= 1
            while arr0_tmp[i_hi_bound] >= LOW_DENSITY_THRESHOLD:
                i_hi_bound += 1
            i_center = int( 0.5 * (i_lo_bound + i_hi_bound) )
            new_index = [k for k in range(i_center - 50, i_center + 50)]
            data0_shift = arr0_tmp[new_index]
            data1_shift = arr1_tmp[new_index]
            pro0_new_data[j, :] = data0_shift[:]
            pro1_new_data[j, :] = data1_shift[:]

            # calculate high density and low density
            lo_dens_region = [k for k in range(i_center - 50, i_lo_bound - 2)] + [k for k in range(i_hi_bound + 3, i_center + 50)]
            hi_dens_region = [k for k in range(i_lo_bound - 2, i_hi_bound + 3)]
            pro0_lo_density = np.mean(arr0_tmp[lo_dens_region])
            pro0_hi_density = np.mean(arr0_tmp[hi_dens_region])
            pro1_lo_density = np.mean(arr1_tmp[lo_dens_region])
            pro1_hi_density = np.mean(arr1_tmp[hi_dens_region])

            # output to density profile
            for k in range(num_bins):
                pro0_dpf.write("{0:18.12f}  ".format( data0_shift[k] ))
                pro1_dpf.write("{0:18.12f}  ".format( data1_shift[k] ))
            pro0_dpf.write("\n")
            pro1_dpf.write("\n")
            # output to density average
            pro0_daf.write("{0:8d}   {1:8d}   {2:8d}   {3:18.12f}   {4:18.12f}  \n".format(i_center, i_lo_bound, i_hi_bound, pro0_lo_density, pro0_hi_density))
            pro1_daf.write("{0:8d}   {1:8d}   {2:8d}   {3:18.12f}   {4:18.12f}  \n".format(i_center, i_lo_bound, i_hi_bound, pro1_lo_density, pro1_hi_density))

        pro0_dpf.close()
        pro0_daf.close()
        pro1_dpf.close()
        pro1_daf.close()

        axes[i * 2].imshow(pro0_new_data.T, cmap="Blues", vmin=0.000001, vmax=.005, origin="lower")
        axes[i * 2 + 1].imshow(pro1_new_data.T, cmap="Reds", vmin=0.000001, vmax=.005, origin="lower")

        axes[i * 2].set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
        axes[i * 2].set_xlim(0, 500)
        axes[i * 2 + 1].set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
        axes[i * 2 + 1].set_xlim(0, 500)

        axes[i * 2].set_yticks([0, 50, 100])
        axes[i * 2].set_yticklabels([0, 15, 30], fontsize=10)
        axes[i * 2].set_ylim(0, 100)
        axes[i * 2].set_ylabel(r"z ($10^2\AA$)", fontsize=12)
        axes[i * 2 + 1].set_yticks([0, 50, 100])
        axes[i * 2 + 1].set_yticklabels([0, 15, 30], fontsize=10)
        axes[i * 2 + 1].set_ylim(0, 100)
        axes[i * 2 + 1].set_ylabel(r"z ($10^2\AA$)", fontsize=12)

        sub_leg_0 = "T={0:>3d} K, {1}".format(t, pro0_name)
        sub_leg_1 = "T={0:>3d} K, {1}".format(t, pro1_name)
        axes[i * 2].text(100, 75, sub_leg_0)
        axes[i * 2 + 1].text(100, 75, sub_leg_1)

    axes[num_temperatures * 2 - 1].set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axes[num_temperatures * 2 - 1].set_xticklabels([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], fontsize=10)
    axes[num_temperatures * 2 - 1].set_xlim(0, 500)
    axes[num_temperatures * 2 - 1].set_xlabel(r"MD steps ($\times 10^6$)", fontsize=12)

    figname = "{0}_{1}_md4_T_density_centered.png".format(pro0_name, pro1_name)
    plt.savefig(figname, dpi=150)
    # plt.show()

if __name__ == '__main__':
    main()
