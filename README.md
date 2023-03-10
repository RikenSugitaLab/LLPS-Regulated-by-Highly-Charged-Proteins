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

We put all the MD simulation files and data analysis scripts in the `simulation_files` directory.  With these files, one can reproduce all the results reported in the paper. MD trajectories are not included.

- `simulation_files/`
  - `01_all_atom_simulation_of_Hero11/`
    - `README.md`
    - `01_WT_from_helical_structure/`
    - `02_WT_from_random_coil/`
    - `03_KRless_from_helical_structure/`
    - `04_KRless_from_random_coil/`
    - `ANALYSIS`
  - `02_CG_single_chain_TDP43LCD`
    - `01_with_alpha`
    - `01_without_alpha`
  - `03_CG_single_chain_Hero11WT`
    - `01_with_alpha`
    - `01_without_alpha`
  - `04_CG_multi_chain_TDP43LCD_homotypic`
    - `01_TDP43LCD_with_alpha`
    - `02_TDP43LCD_without_alpha`
  - `05_CG_multi_chain_Hero11WT_homotypic`
    - `01_Hero11WT_with_alpha`
    - `02_Hero11WT_without_alpha`
  - `06_CG_multi_chain_Hero11_mutants_homotypic`
    - `01.1_KRless_with_alpha`
    - `01.2_KRless_without_alpha`
    - `02.1_chargeless_with_alpha`
    - `02.2_chargeless_without_alpha`
    - `03_scramble_without_alpha`
  - `07_CG_multi_chain_TDP43LCD_100_Hero11WT_100_heterotypic`
    - `01_Hero11WT_with_alpha`
    - `02_Hero11WT_without_alpha`
  - `08_CG_multi_chain_TDP43LCD_100_Hero11mutants_100_heterotypic`
    - `01.1_KRless_with_helix`
    - `01.2_KRless_without_helix`
    - `02.1_chargeless_with_helix`
    - `02.2_chargeless_without_helix`
  - `09_CG_multi_chain_TDP43LCD_100_Hero11WT_N_heterotypic`
    - `01_Hero11WT_with_alpha_multi_temperature`
    - `02_Hero11WT_without_alpha_multi_temperature`
    - `03_Hero11WT_290K_more_simulations`
