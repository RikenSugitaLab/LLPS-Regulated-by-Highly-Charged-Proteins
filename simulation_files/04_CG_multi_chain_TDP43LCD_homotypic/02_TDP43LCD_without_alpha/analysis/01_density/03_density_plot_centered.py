#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():
    T_list = [260 + i * 10 for i in range(10)]
    num_temperatures = len(T_list)
    num_run = 1

    num_steps = 5000
    step_interval = 5
    plot_steps = num_steps // step_interval

    LOW_DENSITY_THRESHOLD = 0.003
    num_bins = 100


    fig, axes = plt.subplots(num_temperatures, 1, figsize=(9, num_temperatures), constrained_layout=True, sharex=True, sharey=False)
    for i, t in enumerate(T_list):
        print("Temperature:", t)
        fname = "tdp43_md2_T_{0:0>2d}.density.dat".format(i + 1)
        data_local = np.loadtxt(fname)[::step_interval, :]

        density_profile_output_fname = "tdp43_md2_T_{0:0>2d}_density_centered.dat".format(i + 1)
        density_average_output_fname = "tdp43_md2_T_{0:0>2d}_density_averaged.dat".format(i + 1)
        dpf = open(density_profile_output_fname, "w")
        daf = open(density_average_output_fname, "w")

        # ================
        # guess the center
        # ================
        new_data = np.zeros((plot_steps, num_bins))
        for j in range(plot_steps):
            # find the peak
            data_j = data_local[j, :]
            arr_tmp = np.append(data_j, data_j)
            arr_tmp = np.append(data_j, arr_tmp)
            i_peak = np.argmax(data_j) + num_bins
            i_lo_bound, i_hi_bound = i_peak, i_peak
            while arr_tmp[i_lo_bound] >= LOW_DENSITY_THRESHOLD:
                i_lo_bound -= 1
            while arr_tmp[i_hi_bound] >= LOW_DENSITY_THRESHOLD:
                i_hi_bound += 1
            i_center = int( 0.5 * (i_lo_bound + i_hi_bound) )
            new_index = [k for k in range(i_center - 50, i_center + 50)]
            data_shift = arr_tmp[new_index]
            new_data[j, :] = data_shift[:]

            # calculate high density and low density
            lo_dens_region = [k for k in range(i_center - 50, i_lo_bound - 2)] + [k for k in range(i_hi_bound + 3, i_center + 50)]
            hi_dens_region = [k for k in range(i_lo_bound - 2, i_hi_bound + 3)]
            lo_density = np.mean(arr_tmp[lo_dens_region])
            hi_density = np.mean(arr_tmp[hi_dens_region])

            # output to density profile
            for k in range(num_bins):
                dpf.write("{0:18.12f}  ".format( data_shift[k] ))
            dpf.write("\n")
            # output to density average
            daf.write("{0:8d}   {1:8d}   {2:8d}   {3:18.12f}   {4:18.12f}  \n".format(i_center, i_lo_bound, i_hi_bound, lo_density, hi_density))

        dpf.close()
        daf.close()

        axes[i].imshow(new_data.T, cmap="Blues", vmin=0.000001, vmax=.005, origin="lower")

        axes[i].set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
        axes[i].set_xlim(0, 1000)

        axes[i].set_yticks([0, 50, 100])
        axes[i].set_yticklabels([0, 15, 30], fontsize=10)
        axes[i].set_ylim(0, 100)
        axes[i].set_ylabel(r"z ($10^2\AA$)", fontsize=12)

        sub_leg = "T={0:>3d} K".format(t)
        axes[i].text(100, 75, sub_leg)

    axes[num_temperatures - 1].set_xticks([0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    axes[num_temperatures - 1].set_xticklabels([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50], fontsize=10)
    axes[num_temperatures - 1].set_xlim(0, 1000)
    axes[num_temperatures - 1].set_xlabel(r"MD steps ($\times 10^6$)", fontsize=12)

    figname = "tdp43_md2_T_density_centered.png".format(t)
    plt.savefig(figname, dpi=150)
    # plt.show()

if __name__ == '__main__':
    main()
