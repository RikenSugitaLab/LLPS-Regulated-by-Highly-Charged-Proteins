mol new ../../01_build/ionize.psf

set hero [atomselect top "protein"]
set idxList [$hero get index]

set idxFile [open hero_protein.ind w]
puts $idxFile $idxList

exit
