mol new ../../01_build/ionize.psf
mol addfile ../../01_build/ionize.pdb

set hero [atomselect top "protein"]
$hero writepdb hero_0.pdb

exit
