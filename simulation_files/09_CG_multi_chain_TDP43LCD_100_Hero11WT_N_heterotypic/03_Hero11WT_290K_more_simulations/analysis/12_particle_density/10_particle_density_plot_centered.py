#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def compute_plot_density(conc_hero11, i_run):
    pro0_name = "tdp43"
    pro1_name = "hero11"

    # T_list = [260 + i * 10 for i in range(7)]
    T_list = [290]
    num_temperatures = len(T_list)
    num_run = 1

    num_steps = 3000
    step_interval = 1
    # step_interval_plot = 3
    # plot_steps = num_steps // step_interval_plot

    LOW_DENSITY_THRESHOLD = 0.003
    num_bins = 100


    # fig, axes = plt.subplots(num_temperatures, num_run * 1, figsize=(num_run * 9, num_temperatures), constrained_layout=True, sharex=True, sharey=False)
    for i, t in enumerate(T_list):
        system_name = "{0}_{1}_{2}_{3}".format(pro0_name, 100, pro1_name, conc_hero11)
        print("Temperature:", t)

        md3_fname0 = "../01_density_in_box/md3/{0}_T{1:0>3d}_md3_r{2:0>2d}.{3}.density.dat".format(system_name, t, i_run, pro0_name)
        md3_fname1 = "../01_density_in_box/md3/{0}_T{1:0>3d}_md3_r{2:0>2d}.{3}.density.dat".format(system_name, t, i_run, pro1_name)
        md3_data0_local = np.loadtxt(md3_fname0)[::step_interval, :]
        md3_data1_local = np.loadtxt(md3_fname1)[::step_interval, :]

        # data0_all = md1_data0_local
        # data1_all = md1_data1_local
        # data0_all = np.vstack((md1_data0_local, md3_data0_local))
        # data1_all = np.vstack((md1_data1_local, md3_data1_local))
        # data0_all = np.vstack((md1_data0_local, md3_data0_local, md3_data0_local))
        # data1_all = np.vstack((md1_data1_local, md3_data1_local, md3_data1_local))
        data0_all = md3_data0_local
        data1_all = md3_data1_local

        pro_all_density_profile_output_fname = "{0}_T_{1:0>3d}_{2}_{3:0>2d}_density_centered.dat".format(system_name, t, "all", i_run)
        pro_all_density_average_output_fname = "{0}_T_{1:0>3d}_{2}_{3:0>2d}_density_averaged.dat".format(system_name, t, "all", i_run)
        # pro_all_dpf = open(pro_all_density_profile_output_fname, "w")
        pro_all_daf = open(pro_all_density_average_output_fname, "w")

        # ================
        # guess the center
        # ================
        # pro0_plot_data = np.zeros((plot_steps, num_bins))
        # pro1_plot_data = np.zeros((plot_steps, num_bins))
        # pro_all_plot_data = np.zeros((plot_steps, num_bins))
        for j in range(num_steps):
            if j >= data0_all.shape[0]:
                break
            # find the peak
            data0_j = data0_all[j, :]
            data1_j = data1_all[j, :]
            arr0_tmp = np.append(data0_j, data0_j)
            arr0_tmp = np.append(data0_j, arr0_tmp)
            arr1_tmp = np.append(data1_j, data1_j)
            arr1_tmp = np.append(data1_j, arr1_tmp)
            i_peak = np.argmax(data0_j) + num_bins
            i_lo_bound, i_hi_bound = i_peak, i_peak
            while arr0_tmp[i_lo_bound] + arr1_tmp[i_lo_bound] >= LOW_DENSITY_THRESHOLD:
                i_lo_bound -= 1
            while arr0_tmp[i_hi_bound] + arr1_tmp[i_hi_bound] >= LOW_DENSITY_THRESHOLD:
                i_hi_bound += 1
            i_center = int( 0.5 * (i_lo_bound + i_hi_bound) )
            new_index = [k for k in range(i_center - 50, i_center + 50)]
            data0_shift = arr0_tmp[new_index]
            data1_shift = arr1_tmp[new_index]
            data_all_shift = data0_shift * 154 + data1_shift * 99
            # if j % step_interval_plot == 1:
                # pro0_plot_data[j // step_interval_plot, :] = data0_shift[:]
                # pro1_plot_data[j // step_interval_plot, :] = data1_shift[:]
                # pro_all_plot_data[j // step_interval_plot, :] = data_all_shift[:]

            # calculate high density and low density
            lo_dens_region = [k for k in range(i_center - 50, i_lo_bound - 2)] + [k for k in range(i_hi_bound + 3, i_center + 50)]
            hi_dens_region = [k for k in range(i_lo_bound - 2, i_hi_bound + 3)]
            # pro0_lo_density = np.mean(arr0_tmp[lo_dens_region])
            # pro0_hi_density = np.mean(arr0_tmp[hi_dens_region])
            # pro1_lo_density = np.mean(arr1_tmp[lo_dens_region])
            # pro1_hi_density = np.mean(arr1_tmp[hi_dens_region])
            pro_all_lo_density = np.mean(arr0_tmp[lo_dens_region] * 154 + arr1_tmp[lo_dens_region] * 99)
            pro_all_hi_density = np.mean(arr0_tmp[hi_dens_region] * 154 + arr1_tmp[hi_dens_region] * 99)


            # output to density profile
            # for k in range(num_bins):
                # pro0_dpf.write("{0:18.12f}  ".format( data0_shift[k] ))
                # pro1_dpf.write("{0:18.12f}  ".format( data1_shift[k] ))
                # pro_all_dpf.write("{0:18.12f}  ".format( data_all_shift[k] ))
            # pro0_dpf.write("\n")
            # pro1_dpf.write("\n")
            # pro_all_dpf.write("\n")
            # output to density average
            # pro0_daf.write("{0:8d}   {1:8d}   {2:8d}   {3:18.12f}   {4:18.12f}  \n".format(i_center, i_lo_bound, i_hi_bound, pro0_lo_density, pro0_hi_density))
            # pro1_daf.write("{0:8d}   {1:8d}   {2:8d}   {3:18.12f}   {4:18.12f}  \n".format(i_center, i_lo_bound, i_hi_bound, pro1_lo_density, pro1_hi_density))
            pro_all_daf.write("{0:8d}   {1:8d}   {2:8d}   {3:18.12f}   {4:18.12f}  \n".format(i_center, i_lo_bound, i_hi_bound, pro_all_lo_density, pro_all_hi_density))

        # pro0_dpf.close()
        # pro0_daf.close()
        # pro1_dpf.close()
        # pro1_daf.close()
        # pro_all_dpf.close()
        pro_all_daf.close()

        # axes[i * 2].imshow(pro0_plot_data.T, cmap="Blues", vmin=0.000001, vmax=.005, origin="lower")
        # axes[i * 2 + 1].imshow(pro1_plot_data.T, cmap="Reds", vmin=0.000001, vmax=.005, origin="lower")
        # axes[i].imshow(pro_all_plot_data.T, cmap="Greens", vmin=0.000001, vmax=.005, origin="lower")

        # axes[i * 2].set_xticks([500 * j for j in range(11)])
        # axes[i].set_xticks([500 * j for j in range(11)])
        # axes[i].set_xlim(0, 3000)
        # axes[i * 2 + 1].set_xticks([500 * j for j in range(11)])
        # axes[i * 2 + 1].set_xlim(0, 900)
        # axes[i * 2 + 1].set_xlim(0, 3000)

        # axes[i].set_yticks([0, 50, 100])
        # axes[i].set_yticklabels([0, 15, 30], fontsize=10)
        # axes[i].set_ylim(0, 100)
        # axes[i].set_ylabel(r"z ($10^2\AA$)", fontsize=12)
        # axes[i * 2].set_yticks([0, 50, 100])
        # axes[i * 2].set_yticklabels([0, 15, 30], fontsize=10)
        # axes[i * 2].set_ylim(0, 100)
        # axes[i * 2].set_ylabel(r"z ($10^2\AA$)", fontsize=12)
        # axes[i * 2 + 1].set_yticks([0, 50, 100])
        # axes[i * 2 + 1].set_yticklabels([0, 15, 30], fontsize=10)
        # axes[i * 2 + 1].set_ylim(0, 100)
        # axes[i * 2 + 1].set_ylabel(r"z ($10^2\AA$)", fontsize=12)

        # sub_leg_0 = "T={0:>3d} K, {1}".format(t, pro0_name)
        # sub_leg_1 = "T={0:>3d} K, {1}".format(t, pro1_name)
        # axes[i * 2].text(100, 75, sub_leg_0)
        # axes[i * 2 + 1].text(100, 75, sub_leg_1)
        # sub_leg_0 = "T={0:>3d} K, {1}".format(t, pro0_name)
        # axes[i].text(100, 75, sub_leg_0)

    # axes[num_temperatures * 2 - 1].set_xticks([500 * j for j in range(11)])
    # axes[num_temperatures * 2 - 1].set_xticklabels([15 * j for j in range(11)], fontsize=10)
    # axes[num_temperatures * 2 - 1].set_xlim(0, 900)
    # axes[num_temperatures * 2 - 1].set_xlim(0, 3000)
    # axes[num_temperatures * 2 - 1].set_xlabel(r"MD steps ($\times 10^6$)", fontsize=12)
    # axes[num_temperatures - 1].set_xticks([500 * j for j in range(11)])
    # axes[num_temperatures - 1].set_xticklabels([15 * j for j in range(11)], fontsize=10)
    # axes[num_temperatures - 1].set_xlim(0, 3000)
    # axes[num_temperatures - 1].set_xlabel(r"MD steps ($\times 10^6$)", fontsize=12)

    # figname = "{0}_T_density_centered_all_particles.png".format(system_name)
    # plt.savefig(figname, dpi=150)
    # plt.show()

if __name__ == '__main__':
    for j in range(1, 10):
        compute_plot_density(10 * j, 1)
        compute_plot_density(10 * j, 2)
        compute_plot_density(10 * j, 3)
        compute_plot_density(10 * j, 4)
        compute_plot_density(10 * j, 5)