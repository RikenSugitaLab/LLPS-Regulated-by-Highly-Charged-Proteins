#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def calculate_1d_contact_count_from_matrix(n_hero, temperature, i_run, cntct_pro1, cntct_pro2, phase_state):
    """Integrate contact matrix to get 1d contact count.
    - n_hero: number of heros
    - temperature: simulation T
    - i_run: #md
    - cntct_pro1: name of protein 1 ("tdp43" or "hero11")
    - cntct_pro2: name of protein 2 ("tdp43" or "hero11")
    - phase_state: dense or dilute phase?
    """

    name_pro1 = "tdp43"
    name_pro2 = "hero11"

    system_name = "{0}_{1}_{2}_{3}".format(name_pro1, 100, name_pro2, n_hero)

    fname = "./md{0}/{1}_t{2:3d}_md{0}.{3}-{4}.{5}.contact_matrix.dat".format(i_run, system_name, temperature, cntct_pro1, cntct_pro2, phase_state)
    contact_data = np.loadtxt(fname)
    zmax = np.max(contact_data)
    print("Maximum of matrix:", zmax)

    contact_num_1d = np.sum(contact_data, axis=1)
    fout_name = "{0}_1d_contact_count.{0}-{1}.{2}.T{3:3d}.{4}_{5}.dat".format(cntct_pro1, cntct_pro2, phase_state, temperature, name_pro2, n_hero)
    fout = open(fout_name, "w")
    for j in range(len(contact_num_1d)):
        fout.write("{0:>3d}    {1:12.6f} \n".format(j + 1, contact_num_1d[j]))

if __name__ == '__main__':
    for T in [260, 270, 280, 290]:
        for i in range(9):
            n = i * 10 + 10
            calculate_1d_contact_count_from_matrix(n, T, 3, "tdp43",  "tdp43",  "dense")
            calculate_1d_contact_count_from_matrix(n, T, 3, "hero11", "hero11", "dense")
            calculate_1d_contact_count_from_matrix(n, T, 3, "tdp43",  "hero11", "dense")
            calculate_1d_contact_count_from_matrix(n, T, 3, "hero11", "tdp43",  "dense")
            calculate_1d_contact_count_from_matrix(n, T, 3, "tdp43",  "tdp43",  "dilute")
            calculate_1d_contact_count_from_matrix(n, T, 3, "hero11", "hero11", "dilute")
            calculate_1d_contact_count_from_matrix(n, T, 3, "tdp43",  "hero11", "dilute")
            calculate_1d_contact_count_from_matrix(n, T, 3, "hero11", "tdp43",  "dilute")
