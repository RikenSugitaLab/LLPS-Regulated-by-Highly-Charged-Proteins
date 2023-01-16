#!/bin/bash

j=1
k=0

for t in {01..10}; do        # temperature
    for r in {01..05}; do        # temperature
        mytemperature=$(echo "$t * 10 + 200" | bc)
        cp pro.atin md${j}_T_${t}_$r.atin
        cp pro.job_hokusai  md${j}_T_${t}_$r.job
        cp hero11_noalpha_100.top  hero11_T_${t}_$r.top
        mycrdfile=$(echo "$r % 5 + 1" | bc)
        cp crd/hero11_bs_final_0$mycrdfile.gro crd/hero11_T_${t}_$r.gro
        sed -e "s/TOPFILE/hero11_T_${t}_$r/g" -i md${j}_T_${t}_$r.atin
        sed -e "s/CRDFILE/hero11_T_${t}_$r/g" -i md${j}_T_${t}_$r.atin
        sed -e "s/TEMPLATE/hero11_noalpha_100_md${j}_T_${t}_$r/g" -i md${j}_T_${t}_$r.atin
        sed -e "/temperature/s/300/$mytemperature/g" -i md${j}_T_${t}_$r.atin
        sed -e "s/SEEDNUM/${j}${r}$t/g" -i md${j}_T_${t}_$r.atin
        sed -e "s/RUNNAME/md${j}_T_${t}_$r/g" -i md${j}_T_${t}_$r.job
        sleep 0.1
        pjsub md${j}_T_${t}_$r.job
    done
done
