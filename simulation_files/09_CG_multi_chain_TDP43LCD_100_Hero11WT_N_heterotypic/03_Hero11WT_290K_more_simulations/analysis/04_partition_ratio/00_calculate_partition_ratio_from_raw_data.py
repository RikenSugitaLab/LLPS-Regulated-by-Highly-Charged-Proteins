#!/usr/bin/env python3

import numpy as np

def calculate_partition_from_density(temperature, i_conc):
    num_runs = 5
    num_md   = 3

    pro1_name = "tdp43"
    pro2_name = "hero11"

    system_name = "{0}_{1}_{2}_{3}".format(pro1_name, 100, pro2_name, i_conc * 10)
    print("analyzing system: ", system_name)

    pro1_density_threshold = 0.0034

    for j in range(num_runs):
        print("   >> run ", j + 1)
        pro2_partition_out = open("{0}_T{1:>03d}_r{2:>02d}.{3}.partition.dat".format(system_name, temperature, j + 1, pro2_name), "w")
        i_step = 0
        for k in range(num_md):
            pro1_fname = "../01_density_in_box/md{2}/{0}_T{1:>03d}_md{2}_r{3:>02d}.{4}.density.dat".format(system_name, temperature, k + 1, j+1, pro1_name)
            pro2_fname = "../01_density_in_box/md{2}/{0}_T{1:>03d}_md{2}_r{3:>02d}.{4}.density.dat".format(system_name, temperature, k + 1, j+1, pro2_name)
            pro1_data = np.loadtxt(pro1_fname)
            pro2_data = np.loadtxt(pro2_fname)

            for tstep in range(pro1_data.shape[0]):
                i_step += 1
                p1_data = pro1_data[tstep]
                p2_data = pro2_data[tstep]
                data_size = p1_data.size
                p2_in_count  = sum(p2_data[j] for j in range(data_size) if p1_data[j] >  pro1_density_threshold)
                p2_out_count = sum(p2_data[j] for j in range(data_size) if p1_data[j] <= pro1_density_threshold)
                pro2_partition_out.write("{0:6d}   {1:>12.9f}   {2:>12.9f} \n".format(i_step, p2_in_count, p2_out_count))
        pro2_partition_out.close()

if __name__ == '__main__':
    temperature_list = [290 + i * 5 for i in range(2)]
    num_conc = 9
    for it, temperature in enumerate(temperature_list):
        for i in range(num_conc):
            calculate_partition_from_density(temperature, i + 1)
