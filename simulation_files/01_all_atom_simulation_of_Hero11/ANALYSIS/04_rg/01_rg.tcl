# mol new ../01_extract_protein_only_dcd/hero_0.pdb
mol new ./hero_native.pdb
mol addfile ../01_extract_protein_only_dcd/hero11_WT_helix_md_all.dcd waitfor all

set allca [atomselect 0 "name CA"]
set rgFile [open hero_rg.dat w]

set nf [molinfo 0 get numframes]
for { set i 0 } { $i <= $nf } { inc i } {
    $allca frame $i
    set rgCA [measure rgyr $allca]
    puts $rgFile $rgCA
}

exit

