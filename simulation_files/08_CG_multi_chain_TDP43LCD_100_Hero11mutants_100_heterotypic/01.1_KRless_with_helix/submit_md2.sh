#!/bin/bash

j=2
k=1

for t in {01..10}; do        # epsilon
    mytemperature=$(echo "240 + $t * 10" | bc)
    cp pro_md${j}.atin md${j}_T_$t.atin
    cp pro.job_hokusai  md${j}_T_$t.job
    # cp hero11_100_tdp43_ctf_100.top  ht_T_$t.top
    # mycrdfile=$(echo "$t % 5 + 1" | bc)
    # cp crd/md1.gro crd/ht_T_$t.gro
    sed -e "s/TOPFILE/ht_T_$t/g" -i md${j}_T_$t.atin
    sed -e "s/CRDFILE/ht_T_$t/g" -i md${j}_T_$t.atin
    sed -e "s/RSTFILE/hero11_tdp43_md${k}_T_$t/g" -i md${j}_T_$t.atin
    sed -e "s/TEMPLATE/hero11_tdp43_md${j}_T_$t/g" -i md${j}_T_$t.atin
    sed -e "/temperature/s/300/$mytemperature/g" -i md${j}_T_$t.atin
    sed -e "s/SEEDNUM/${j}5$t/g" -i md${j}_T_$t.atin
    sed -e "s/RUNNAME/md${j}_T_$t/g" -i md${j}_T_$t.job
    sleep 0.1
    pjsub md${j}_T_$t.job
done

