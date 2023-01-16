#!/bin/bash

j=1
k=0

for t in {01..20}; do        # epsilon
    mytemperature=$(echo "$t * 10 + 150" | bc)
    cp pro.atin md${j}_T_$t.atin
    cp pro.job_eel04  md${j}_T_$t.job
    cp hero11_100_tdp43_ctf_100.top  ht_T_$t.top
    mycrdfile=$(echo "$t % 5 + 1" | bc)
    cp crd/hero11_tdp43_bs_final_0$mycrdfile.gro crd/ht_T_$t.gro
    sed -e "s/TOPFILE/ht_T_$t/g" -i md${j}_T_$t.atin
    sed -e "s/CRDFILE/ht_T_$t/g" -i md${j}_T_$t.atin
    sed -e "s/TEMPLATE/hero11_tdp43_md${j}_T_$t/g" -i md${j}_T_$t.atin
    sed -e "/temperature/s/300/$mytemperature/g" -i md${j}_T_$t.atin
    sed -e "s/SEEDNUM/${j}1$t/g" -i md${j}_T_$t.atin
    sed -e "s/RUNNAME/md${j}_T_$t/g" -i md${j}_T_$t.job
    sleep 0.1
    qsub md${j}_T_$t.job
done

