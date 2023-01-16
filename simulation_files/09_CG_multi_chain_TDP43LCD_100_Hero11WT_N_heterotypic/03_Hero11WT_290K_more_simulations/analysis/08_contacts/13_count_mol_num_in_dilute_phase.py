#!/usr/bin/env python3

import numpy as np

def count_dilute_contacts(T, n_hero):
    n_run = 5

    cnt_00_all = []
    cnt_11_all = []
    cnt_sum_all = []
    for i_run in range(n_run):
        # hero11-hero11
        cnt_00 = []
        fname0 = "./md2/tdp43_100_hero11_{2}_T{0}_md2_r{1:0>2d}.hero11-hero11.dilute.contact_count_ts.dat"
        fname = fname0.format(T, i_run + 1, n_hero)
        with open(fname, "r") as fin:
            for line in fin:
                tmp_sum = 0
                words = line.split()
                if len(words) < 1:
                    continue
                for w in words:
                    cnt = float(w)
                    if cnt >= -0.5:
                        tmp_sum += 1
                cnt_00.append(tmp_sum)

        # tdp43-tdp43
        cnt_11 = []
        fname0 = "./md2/tdp43_100_hero11_{2}_T{0}_md2_r{1:0>2d}.tdp43-tdp43.dilute.contact_count_ts.dat"
        fname = fname0.format(T, i_run + 1, n_hero)
        with open(fname, "r") as fin:
            for line in fin:
                tmp_sum = 0
                words = line.split()
                if len(words) < 1:
                    continue
                for w in words:
                    cnt = float(w)
                    if cnt >= -0.5:
                        tmp_sum += 1
                cnt_11.append(tmp_sum)

        # sum
        cnt_00_arra = np.array(cnt_00)
        cnt_11_arra = np.array(cnt_11)
        cnt_sum_arra = cnt_00_arra + cnt_11_arra

        cnt_00_all.append(cnt_00_arra)
        cnt_11_all.append(cnt_11_arra)
        cnt_sum_all.append(cnt_sum_arra)

    # calculate mean, std
    cnt_00_flat = np.ravel(cnt_00_all)
    cnt_11_flat = np.ravel(cnt_11_all)
    cnt_sum_flat = np.ravel(cnt_sum_all)

    cnt_00_mean = np.mean(cnt_00_flat)
    cnt_00_std  = np.std(cnt_00_flat)
    cnt_00_err  = cnt_00_std / np.sqrt(len(cnt_00_flat))
    cnt_11_mean = np.mean(cnt_11_flat)
    cnt_11_std  = np.std(cnt_11_flat)
    cnt_11_err  = cnt_11_std / np.sqrt(len(cnt_11_flat))
    cnt_sum_mean = np.mean(cnt_sum_flat)
    cnt_sum_std  = np.std(cnt_sum_flat)
    cnt_sum_err  = cnt_sum_std / np.sqrt(len(cnt_sum_flat))

    fout_name = "MOL_NUM_T{0}_nHero_{1}_dilute_all_MD2_REVISION.dat".format(T, n_hero)
    fout = open(fout_name, "w")
    fout.write("ALL: {0:12.3f}   {1:12.3f}   {2:12.3f} \n".format(cnt_sum_mean, cnt_sum_std, cnt_sum_err))
    fout.write("000: {0:12.3f}   {1:12.3f}   {2:12.3f} \n".format(cnt_00_mean, cnt_00_std, cnt_00_err))
    fout.write("111: {0:12.3f}   {1:12.3f}   {2:12.3f} \n".format(cnt_11_mean, cnt_11_std, cnt_11_err))
    fout.close()

if __name__ == '__main__':
    for n_hero in [i* 10 + 10 for i in range(9)]:
        count_dilute_contacts(290, n_hero)
