#!/usr/bin/env python3

import numpy as np

def calculate_partition_from_density(temperature, i_run):
    num_runs = 3

    pro1_name = "tdp43"
    pro2_name = "hero11"

    system_name = "{0}_{1}_{2}_{3}".format(pro1_name, 100, pro2_name, i_run * 10)
    print("analyzing system: ", system_name)

    pro1_density_threshold = 0.0034
    pro2_partition = []

    pro2_partition_out = open("{1}_t{0:>03d}.{2}.partition.dat".format(temperature, system_name, pro2_name), "w")

    i_step = 0
    for j in range(num_runs):
        pro1_fname = "../01_density_in_box/md{3}/{1}_md{3}_T_{0:>03d}.{2}.density.dat".format(temperature, system_name, pro1_name, j+1)
        pro2_fname = "../01_density_in_box/md{3}/{1}_md{3}_T_{0:>03d}.{2}.density.dat".format(temperature, system_name, pro2_name, j+1)
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
    temperature_list = [260 + i * 10 for i in range(4)]
    num_rep = 9
    for it, temperature in enumerate(temperature_list):
        for i in range(num_rep):
            calculate_partition_from_density(temperature, i + 1)
