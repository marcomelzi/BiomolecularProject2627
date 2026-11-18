# Ricalcola solo il q-value separato per classe, su una traiettoria gia' esistente.
# Lanciare DENTRO una cartella run*, con:
#   vmd -dispdev text -e /Users/User/Desktop/BiomolecularProject2627/run_split.tcl
mol new ionized.psf type psf waitfor all
mol addfile prot3.dcd type dcd waitfor all
source /Users/User/Desktop/BiomolecularProject2627/nativecontacts_split.tcl
exit
