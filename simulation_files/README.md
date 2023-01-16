## All-atom simulations

See details in [01_all_atom_simulation_of_Hero11](https://github.com/RikenSugitaLab/LLPS-Regulated-by-Highly-Charged-Proteins/blob/main/simulation_files/01_all_atom_simulation_of_Hero11/README.md).

## CG simulations

Please copy the `param` subdirectory to each system before you start the simulations.  The `param` directory stores CG parameter files.  It is the same as the one in [GENESIS-CG-tool](https://github.com/genesis-release-r-ccs/genesis_cg_tool/tree/main/param), except for the HPS parameters in file `atom_types.itp`, from which we used the HPS-Urry set in the current study.

Basically, to perform CG simulations with GENESIS `atdyn`, please execute the following commands (just an example):
```bash
export OMP_NUM_THREADS=2
mpirun -np 8 -ppn 8 GENESIS/bin/atdyn xxx.atin > xxx.log
```
where `xxx.atin` is the MD control file you can find in each directory.

In many cases, to easily change simulation conditions such as temperature, I wrote bash scripts to generate files automatically.  These files are usually named as "`submit.sh`". 

### Analysis

In the CG simulation directories, there are `analysis` subdirectories, in which we put all of our trajectory analysis scripts.

Many of my scripts calls [GENESIS-CG-tool](https://github.com/genesis-release-r-ccs/genesis_cg_tool).  Please refer to [the wikipage](https://github.com/genesis-release-r-ccs/genesis_cg_tool/wiki) for details.


