# build_system.tcl  --  da PDB (anche da FoldX) a sistema solvatato e ionizzato per NAMD
# Replica il primo capitolo del tutorial NAMD (psfgen -> solvate) e aggiunge autoionize.
#
# uso:
#   vmd -dispdev text -eofexit -e build_system.tcl -args <in.pdb> <outdir> <topologia.inp> <padding_A> <conc_M> <segname>
#
# conc_M = 0 -> nessuna ionizzazione (solo solvatazione, come nel tutorial)

package require psfgen
package require solvate
package require autoionize

lassign $argv in outdir topo pad conc seg
set in   [file normalize $in]
set topo [file normalize $topo]
file mkdir $outdir
cd $outdir

# 1) solo proteina, senza idrogeni (come nel tutorial: atomselect "protein" + writepdb)
#    gli H di FoldX vengono scartati e ricostruiti da psfgen/guesscoord
mol new $in
set sel [atomselect top "protein and noh"]
$sel writepdb prot.pdb
$sel delete
mol delete top

# 2) psfgen (identico a ubq.pgn del tutorial)
resetpsf
topology $topo
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD
segment $seg {pdb prot.pdb}
coordpdb prot.pdb $seg
guesscoord
writepdb prot_psf.pdb
writepsf prot_psf.psf

# 3) solvatazione (water box, -t = padding)
solvate prot_psf.psf prot_psf.pdb -t $pad -o solv
file delete -force combine.pdb combine.psf

# 4) ionizzazione: neutralizza e porta alla concentrazione richiesta
if {$conc > 0} {
    autoionize -psf solv.psf -pdb solv.pdb -neutralize -cc $conc -o ionized
    set final ionized
} else {
    set final solv
}

# 5) parametri della cella periodica per il .conf di NAMD (come nel tutorial 1.5)
mol new $final.psf
mol addfile $final.pdb
set all [atomselect top all]
set mm [measure minmax $all]
lassign $mm lo hi
set f [open box.txt w]
puts $f "cellBasisVector1 [format %.3f [expr {[lindex $hi 0]-[lindex $lo 0]}]] 0 0"
puts $f "cellBasisVector2 0 [format %.3f [expr {[lindex $hi 1]-[lindex $lo 1]}]] 0"
puts $f "cellBasisVector3 0 0 [format %.3f [expr {[lindex $hi 2]-[lindex $lo 2]}]]"
puts $f "cellOrigin [measure center $all]"
close $f
set q [measure sumweights $all weight charge]
puts "== $final: [llength [$all get index]] atomi, carica netta = [format %.3f $q]"
mol delete top
exit
