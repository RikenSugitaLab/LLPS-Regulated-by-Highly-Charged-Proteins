#!/bin/bash

j=3
k=2

for i in {01..09}; do        # epsilon
for t in {01..10}; do        # epsilon
    mynum=$(echo "$i * 10" | bc)
    heronumber=$(printf %03d $mynum)
    mytmp=$(echo "$t * 10 + 250" | bc)
    cp pro.atin md${j}_N_${mynum}_T_$mytmp.atin
    cp pro.job_hokusai  md${j}_${mynum}_$mytmp.job
    cp demix_tdp43_ctf_100_hero11_${heronumber}.top md${j}_N_${mynum}_T_$mytmp.top
    cp crd/demix_tdp43_ctf_100_hero11_${heronumber}.gro crd/md${j}_N_${mynum}_T_$mytmp.gro
    sed -e "s/TOPFILE/md${j}_N_${mynum}_T_$mytmp/g" -i md${j}_N_${mynum}_T_$mytmp.atin
    sed -e "s/CRDFILE/md${j}_N_${mynum}_T_$mytmp/g" -i md${j}_N_${mynum}_T_$mytmp.atin
    sed -e "s/RSTFILE/tdp43_100_hero11_${mynum}_T${mytmp}_md$k/g" -i md${j}_N_${mynum}_T_$mytmp.atin
    sed -e "s/TEMPLATE/tdp43_100_hero11_${mynum}_T${mytmp}_md${j}/g" -i md${j}_N_${mynum}_T_$mytmp.atin
    sed -e "/temperature/s/250/$mytmp/g" -i md${j}_N_${mynum}_T_$mytmp.atin
    sed -e "s/SEEDNUM/${j}${i}$t/g" -i md${j}_N_${mynum}_T_$mytmp.atin
    sed -e "s/RUNNAME/md${j}_N_${mynum}_T_$mytmp/g" -i md${j}_${mynum}_$mytmp.job
    sleep 0.1
    pjsub md${j}_${mynum}_$mytmp.job
done
done

