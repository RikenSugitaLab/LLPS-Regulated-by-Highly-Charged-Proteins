#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def main():
    T_list = [120 + i * 10 for i in range(0, 12)]
    num_run = 5
    num_steps = 500

    num_atom_hero9 = 83
    num_hero9 = 500
    V = 15 * 15 * 2

    for t in T_list:
        print("Temperature:", t)
        fig, axes = plt.subplots(num_run, 1, figsize=(9, num_run), constrained_layout=True, sharex=True, sharey=False)
        for i in range(num_run):
            fname = "hero9_md1_t{0:3d}_{1:0>2d}.density.dat".format(t, i + 1)
            data_local = np.loadtxt(fname) * num_atom_hero9 * num_hero9 / V
            new_index = [j for j in range(62, 100)] + [j for j in range(62)]
            data_local = data_local[:, new_index]
            print(np.max(data_local))

            axes[i].imshow(data_local.T, cmap="Reds", vmin=0.001, vmax=5.00, origin="lower")

            axes[i].set_xticks([0, 200, 400, 600, 800, 1000])
            axes[i].set_xlim(0, 1000)

            axes[i].set_yticks([0, 50, 100])
            axes[i].set_yticklabels([0, 10, 20], fontsize=10)
            axes[i].set_ylim(0, 100)
            axes[i].set_ylabel(r"z ($10^2\AA$)", fontsize=12)

        axes[num_run - 1].set_xticks([0, 200, 400, 600, 800, 1000])
        axes[num_run - 1].set_xticklabels([0, 2, 4, 6, 8, 10], fontsize=10)
        axes[num_run - 1].set_xlim(0, 1000)
        axes[num_run - 1].set_xlabel(r"MD steps ($\times 10^6$)", fontsize=12)

        figname = "hero9_md1_T{0:3d}_densities_evo.svg".format(t)
        plt.savefig(figname, dpi=150)

if __name__ == '__main__':
    main()
