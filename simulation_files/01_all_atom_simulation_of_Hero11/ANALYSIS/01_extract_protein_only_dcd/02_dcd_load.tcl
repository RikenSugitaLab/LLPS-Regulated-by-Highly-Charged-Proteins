mol new ./hero_0.pdb
mol addfile ./hero11_wt_native_md_all.dcd waitfor all

set all [atomselect 0 "all"]
set part [atomselect 0 "resid 40 to 70 and backbone"]
set refconf [atomselect 0 "resid 40 to 70 and backbone" frame 0]

set nf [molinfo 0 get numframes]
for { set i 0 } { $i <= $nf } { inc i } {
    $all frame $i
    $part frame $i
    $all move [measure fit $part $refconf]
}
