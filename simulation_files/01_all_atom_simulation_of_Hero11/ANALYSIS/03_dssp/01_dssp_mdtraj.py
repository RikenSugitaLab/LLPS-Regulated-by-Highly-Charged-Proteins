#!/usr/bin/env python3

import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

def main():
    t = md.load("../01_extract_protein_only_dcd/hero11_WT_helix_md_all.dcd", top="./solute.psf")

    dssp = md.compute_dssp(t)

    dssp_num = np.zeros((dssp.shape[1], dssp.shape[0], 3))

    # =============
    # ss color code
    # =============
    ss_color_code = { "H": (244/255, 64/255, 94/255), "E": (20/255, 184/255, 166/255), "C": (241/255, 245/255, 249/255)}


    # =============
    # output to txt
    # =============
    dssp_of = open("hero11_dssp.dat", "w")
    for i in range(dssp.shape[0]):
        for j in range(dssp.shape[1]):
            dssp_of.write("{0:1}".format(dssp[i, j]))
            s_c = ss_color_code[dssp[i, j]]
            dssp_num[j, i, :] = [s_c[0], s_c[1], s_c[2]]
        dssp_of.write(" \n")
    dssp_of.close()

    # ====
    # plot
    # ====
    fig, ax = plt.subplots(1, 1, figsize=(10, 3), constrained_layout=True, sharex=False, sharey=False)
    ax.imshow(dssp_num, origin="lower", aspect="equal")
    ax.set_xticks([95.238 * i for i in range(11)])
    ax.set_xticklabels([round(0.1 * i, 1) for i in range(11)], fontsize=12)
    ax.set_xlim(0, 952.4)
    ax.set_xlabel(r"time ($\mu s$)", fontsize=16)
    ax.set_yticks([29, 59, 89])
    ax.set_yticklabels([30, 60, 90], fontsize=12)
    ax.set_ylim(0, 98)
    ax.set_ylabel("#res", fontsize=16)

    plt.savefig("hero11_dssp.svg", dpi=300)

if __name__ == '__main__':
    main()
