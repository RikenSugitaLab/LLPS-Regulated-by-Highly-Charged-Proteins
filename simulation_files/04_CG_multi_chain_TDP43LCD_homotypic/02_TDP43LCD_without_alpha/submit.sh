#!/bin/bash

j=1
k=0

for t in {01..20}; do        # epsilon
    mytemperature=$(echo "$t * 5 + 250" | bc)
    cp pro.atin md${j}_T_$t.atin
    cp pro.job_eel10  md${j}_T_$t.job
    cp tdp43_100.top  tdp43_cg_T_$t.top
    mycrdfile=$(echo "$t % 5 + 1" | bc)
    cp crd/tdp43_bs_final_0$mycrdfile.gro crd/tdp43_cg_T_$t.gro
    sed -e "s/TOPFILE/tdp43_cg_T_$t/g" -i md${j}_T_$t.atin
    sed -e "s/CRDFILE/tdp43_cg_T_$t/g" -i md${j}_T_$t.atin
    sed -e "s/TEMPLATE/tdp43_100_md${j}_T_$t/g" -i md${j}_T_$t.atin
    sed -e "/temperature/s/300/$mytemperature/g" -i md${j}_T_$t.atin
    sed -e "s/SEEDNUM/${j}2$t/g" -i md${j}_T_$t.atin
    sed -e "s/RUNNAME/md${j}_T_$t/g" -i md${j}_T_$t.job
    sleep 0.1
    qsub md${j}_T_$t.job
done
