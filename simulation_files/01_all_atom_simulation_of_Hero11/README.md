## All-atom MD simulations of Hero11

In this subdirectory we have all the input files for simulations of Hero11-WT and Hero11-KRless mutant, starting from AlphaFold predicted structure and random coils, respectively.

To reproduce our results, please use [GENESIS 2.0](https://www.r-ccs.riken.jp/labs/cbrt/genesis-version-2-0/).

The force field we used is CHARMM36m [^1].  Due to the copyright, we would suggest the user to download force fields from [this website](https://www.charmm.org/archive/charmm/resources/charmm-force-fields/) and put directory `toppar` in each subdirectory.  One can read `*.inp` (GENESIS MD control file) for more information.


Generally, to run simulations with GENESIS `spdyn`, please execute the following command (just an example):

```bash
export OMP_NUM_THREADS=2
mpirun -np 80 -ppn 20 GENESIS/bin/spdyn xxx.inp > xxx.log
```


## References

[^1]: Huang, J. et al. CHARMM36m: an improved force field for folded and intrinsically disordered proteins. Nat Methods 14, nmeth.4067 (2016).
  
