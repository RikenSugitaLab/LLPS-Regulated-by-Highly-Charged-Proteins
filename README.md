# Highly Charged Proteins and Their Repulsive Interactions Passively Antagonize Biomolecular Condensation

This repository is used to share core data and scripts for the publication of:
- Cheng Tan, Ai Niitsu, and Yuji Sugita (2022) Highly Charged Proteins and Their Repulsive Interactions Passively Antagonize Biomolecular Condensation. *bioRxiv* 10.1101/2022.11.16.516834


## Softwares

### MD simulations

In this work we performed multiscale molecular dynamics (MD) simulations using software GENESIS. It can be downloaded from its website and Github:
- [GENESIS website](https://www.r-ccs.riken.jp/labs/cbrt/genesis-version-2-0/)
- [GENESIS Github repository](https://github.com/genesis-release-r-ccs/genesis)

### Data analysis

For data analysis, we wrote scripts in Julia, Python and Tcl/tk scripts with the help of the following tools:
- CG trajectory analysis: [GENESIS-CG-tool](https://github.com/genesis-release-r-ccs/genesis_cg_tool)
- All-atom trajectory analysis and structural visualization: [VMD](https://www.ks.uiuc.edu/Research/vmd/)
- Protein secondary structure feature analysis: [MDtraj](https://www.mdtraj.org/1.9.8.dev0/index.html)


## Directory Structure

We put all the MD simulation files in the `simulation_files` directory.  With these files, one can reproduce all the results reported in the paper. MD trajectories are not included.

- `simulation_files/`
  - `01_all_atom_simulation_of_Hero11/`
    - `README.md`
    - `01_WT_from_helical_structure/`
    - `02_WT_from_random_coil/`
    - `03_KRless_from_helical_structure/`
    - `04_KRless_from_random_coil/`
