# q_value_per_residue.tcl -- Q per residuo (frazione di contatti nativi conservati)
#
# Stessa definizione di contatto di nativecontacts.tcl / nativecontacts_split.tcl:
#   - contatto = coppia di CA entro 7.0 A, con separazione in sequenza |resid_i - resid_j| >= 3
#   - stato nativo = frame 0 della traiettoria
# Per ogni residuo r:  Q_r(frame) = 100 * (contatti nativi di r ancora presenti) / (contatti nativi di r)
#
# Uso, da dentro una cartella run*  (con ionized.psf e prot3.dcd):
#   vmd -dispdev text -eofexit -e q_value_per_residue.tcl
# Oppure con la molecola gia' caricata (via source da un altro script).
#
# Finestra di media: dal frame 200 all'ultimo (come nella tabella 510_per_residuo.txt).
# Per cambiarla:  QRES_FIRST=100 vmd -dispdev text -eofexit -e q_value_per_residue.tcl
#
# Scrive:
#   q_per_residue.txt         resid resname n_nativi q_res   (q_res = media sulla finestra, in %)
#   q_per_residue_frames.txt  frame + Q (%) di ogni residuo, per ogni frame (nan = nessun contatto nativo)

set pr_cutoff 7.0
set pr_minsep 3
set pr_first  200
if {[info exists env(QRES_FIRST)]} { set pr_first $env(QRES_FIRST) }

# se non c'e' nessuna molecola caricata, la carico io (modalita' standalone)
set pr_standalone 0
if {[molinfo num] == 0} {
    mol new ionized.psf type psf waitfor all
    mol addfile prot3.dcd type dcd waitfor all
    set pr_standalone 1
}

set pr_nf [molinfo top get numframes]
if {$pr_nf == 0} {
    puts "ERRORE: nessun frame letto (psf/dcd mancanti o errati)"
    if {$pr_standalone} { exit } else { error "q_value_per_residue: nessun frame" }
}
if {$pr_first >= $pr_nf} {
    puts "ATTENZIONE: finestra da frame $pr_first ma ci sono solo $pr_nf frame: uso tutti i frame"
    set pr_first 0
}

set pr_sel  [atomselect top "protein and name CA"]
set pr_idx  [$pr_sel get index]
set pr_res  [$pr_sel get resid]
set pr_name [$pr_sel get resname]
if {[llength [lsort -unique $pr_res]] != [llength $pr_res]} {
    puts "ERRORE: resid non univoci tra i CA (piu' catene o codici di inserzione?)"
    if {$pr_standalone} { exit } else { error "q_value_per_residue: resid non univoci" }
}

# mappa atom index -> resid
foreach i $pr_idx r $pr_res { set pr_resid($i) $r }

proc pr_contactSet {sel cutoff minsep} {
    global pr_resid
    set c [measure contacts $cutoff $sel]
    set out {}
    foreach i [lindex $c 0] j [lindex $c 1] {
        if {abs($pr_resid($i) - $pr_resid($j)) < $minsep} continue
        if {$i < $j} { lappend out "$i $j" } else { lappend out "$j $i" }
    }
    return [lsort -unique $out]
}

# --- stato nativo: frame 0 ---
$pr_sel frame 0
set pr_ref [pr_contactSet $pr_sel $pr_cutoff $pr_minsep]
foreach r $pr_res { set pr_nnat($r) 0 ; set pr_sum($r) 0.0 }
foreach p $pr_ref {
    lassign $p i j
    set pr_isnat($p) 1
    incr pr_nnat($pr_resid($i))
    incr pr_nnat($pr_resid($j))
}
puts "Q per residuo: [llength $pr_ref] contatti nativi su [llength $pr_res] residui"

# --- tutti i frame ---
set pr_out [open ./q_per_residue_frames.txt w]
puts $pr_out "# frame, poi Q (%) di ogni residuo nell'ordine indicato sotto (nan = nessun contatto nativo)"
puts $pr_out "# resid $pr_res"
puts $pr_out "# resname $pr_name"

set pr_nwin 0
for {set f 0} {$f < $pr_nf} {incr f} {
    $pr_sel frame $f
    foreach r $pr_res { set pr_nc($r) 0 }
    foreach p [pr_contactSet $pr_sel $pr_cutoff $pr_minsep] {
        if {[info exists pr_isnat($p)]} {
            lassign $p i j
            incr pr_nc($pr_resid($i))
            incr pr_nc($pr_resid($j))
        }
    }
    set riga $f
    foreach r $pr_res {
        if {$pr_nnat($r) > 0} {
            set q [expr {100.0 * $pr_nc($r) / $pr_nnat($r)}]
            if {$f >= $pr_first} { set pr_sum($r) [expr {$pr_sum($r) + $q}] }
            append riga [format " %.3f" $q]
        } else {
            append riga " nan"
        }
    }
    if {$f >= $pr_first} { incr pr_nwin }
    puts $pr_out $riga
    if {$f % 100 == 0} { puts "frame $f" }
}
close $pr_out

# --- tabella riassuntiva (stesse prime 4 colonne di 510_per_residuo.txt) ---
set pr_out [open ./q_per_residue.txt w]
puts $pr_out "# contatti CA-CA entro $pr_cutoff A, separazione >= $pr_minsep, riferimento frame 0"
puts $pr_out "# finestra frame $pr_first-[expr {$pr_nf - 1}] ($pr_nwin frame)"
puts $pr_out "# resid resname n_nativi q_res"
foreach r $pr_res n $pr_name {
    if {$pr_nnat($r) > 0 && $pr_nwin > 0} {
        set q [format "%.2f" [expr {$pr_sum($r) / $pr_nwin}]]
    } else {
        set q nan
    }
    puts $pr_out "$r $n $pr_nnat($r) $q"
}
close $pr_out

if {$pr_standalone} { exit }
