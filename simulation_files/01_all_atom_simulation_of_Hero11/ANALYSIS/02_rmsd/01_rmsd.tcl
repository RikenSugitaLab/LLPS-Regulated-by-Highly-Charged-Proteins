mol new ../01_extract_protein_only_dcd/hero_0.pdb
mol addfile ../01_extract_protein_only_dcd/hero11_WT_helix_md_all.dcd waitfor all

set all [atomselect 0 "resid 38 to 72 and noh"]
set part [atomselect 0 "resid 38 to 72 and backbone"]
set refconf [atomselect 0 "resid 38 to 72 and backbone" frame 0]
set refcomp [atomselect 0 "resid 38 to 72 and noh" frame 0]

set rmsdFile [open hero_rmsd.dat w]

set nf [molinfo 0 get numframes]
for { set i 0 } { $i <= $nf } { inc i } {
    $all frame $i
    $part frame $i
    $all move [measure fit $part $refconf]
    set rmsdVal [measure rmsd $all $refcomp]
    puts $rmsdFile $rmsdVal
}

exit

