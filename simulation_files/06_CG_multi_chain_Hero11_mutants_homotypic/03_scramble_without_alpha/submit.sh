#!/bin/bash

j=1
k=0

for t in {01..08}; do        # epsilon
    for i in {1..3}; do
        mytemperature=$(echo "$t * 10 + 80" | bc)
        cp pro.atin md${j}_scrmb_${i}_T_$t.atin
        cp pro.job_hokusai  md${j}_scrmb_${i}_T_$t.job
        cp hero11_noalpha_100.top  hero11_scrmb_${i}_T_$t.top
        sed -e "s/ITPNAME/hero11_scramble_${i}_noalpha/g" -i hero11_scrmb_${i}_T_$t.top
        mycrdfile=$(echo "$t % 5 + 1" | bc)
        cp crd/hero11_bs_final_0$mycrdfile.gro crd/hero11_scrmb_${i}_T_$t.gro
        sed -e "s/TOPFILE/hero11_scrmb_${i}_T_$t/g" -i md${j}_scrmb_${i}_T_$t.atin
        sed -e "s/CRDFILE/hero11_scrmb_${i}_T_$t/g" -i md${j}_scrmb_${i}_T_$t.atin
        sed -e "s/TEMPLATE/hero11_noalpha_100_md${j}_scrmb_${i}_T_$t/g" -i md${j}_scrmb_${i}_T_$t.atin
        sed -e "/temperature/s/300/$mytemperature/g" -i md${j}_scrmb_${i}_T_$t.atin
        sed -e "s/SEEDNUM/${j}4$t/g" -i md${j}_scrmb_${i}_T_$t.atin
        sed -e "s/RUNNAME/md${j}_scrmb_${i}_T_$t/g" -i md${j}_scrmb_${i}_T_$t.job
        sleep 0.1
        pjsub md${j}_scrmb_${i}_T_$t.job
    done
done
