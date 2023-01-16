#!/bin/bash

j=3
k=2

for i in {01..09}; do        # epsilon
    for t in {01..02}; do        # epsilon
        for irun in {01..05}; do
            mynum=$(echo "$i * 10" | bc)
            heronumber=$(printf %03d $mynum)
            mytmp=$(echo "$t * 5 + 285" | bc)
            cp pro_md${j}.atin md${j}_N_${mynum}_T_${mytmp}_r_${irun}.atin
            cp pro.job_hokusai md${j}_${mynum}_${mytmp}_r_${irun}.job
            cp demix_tdp43_ctf_100_hero11_${heronumber}.top N_${mynum}_T_${mytmp}_r_${irun}.top
            cp crd/demix_tdp43_ctf_100_hero11_${heronumber}.gro crd/N_${mynum}_T_${mytmp}_r_${irun}.gro
            sed -e "s/TOPFILE/N_${mynum}_T_${mytmp}_r_${irun}/g" -i md${j}_N_${mynum}_T_${mytmp}_r_${irun}.atin
            sed -e "s/CRDFILE/N_${mynum}_T_${mytmp}_r_${irun}/g" -i md${j}_N_${mynum}_T_${mytmp}_r_${irun}.atin
            sed -e "s/RSTFILE/tdp43_100_hero11_${mynum}_T${mytmp}_md${k}_r${irun}/g" -i md${j}_N_${mynum}_T_${mytmp}_r_${irun}.atin
            sed -e "s/TEMPLATE/tdp43_100_hero11_${mynum}_T${mytmp}_md${j}_r${irun}/g" -i md${j}_N_${mynum}_T_${mytmp}_r_${irun}.atin
            sed -e "/temperature/s/250/${mytmp}/g" -i md${j}_N_${mynum}_T_${mytmp}_r_${irun}.atin
            sed -e "s/SEEDNUM/${j}${i}${irun}/g" -i md${j}_N_${mynum}_T_${mytmp}_r_${irun}.atin
            sed -e "s/RUNNAME/md${j}_N_${mynum}_T_${mytmp}_r_${irun}/g" -i md${j}_${mynum}_${mytmp}_r_${irun}.job
            sleep 0.1
            pjsub md${j}_${mynum}_${mytmp}_r_${irun}.job
        done
    done
done
