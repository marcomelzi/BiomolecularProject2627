# q-value senza il plugin nativecontacts (solo comandi core di VMD)
set cutoff 7.0
set minsep 3        ;# scarta le coppie con |resid_i - resid_j| < 3

set protein [atomselect top "protein and name CA"]

# mappa atom index -> resid, costruita una volta sola
foreach i [$protein get index] r [$protein get resid] { set resid($i) $r }

proc contactSet {sel cutoff minsep} {
    global resid
    set c [measure contacts $cutoff $sel]
    set out {}
    foreach i [lindex $c 0] j [lindex $c 1] {
        if {abs($resid($i) - $resid($j)) < $minsep} continue
        if {$i < $j} { lappend out "$i $j" } else { lappend out "$j $i" }
    }
    return [lsort -unique $out]
}

set nf [molinfo top get numframes]

## Use first trajectory frame as a reference
$protein frame 0
set ref [contactSet $protein $cutoff $minsep]
set nnc [llength $ref]
foreach p $ref { set isnative($p) 1 }

set out [open ./q_value.txt w]
puts "There are $nnc contacts in the native state"
puts $out "There are $nnc contacts in the native state"

for {set f 0} {$f < $nf} {incr f} {
    $protein frame $f
    set nc 0
    foreach p [contactSet $protein $cutoff $minsep] {
        if {[info exists isnative($p)]} { incr nc }
    }
    set qnc [expr {100.0 * $nc / $nnc}]
    puts       [format "Frame %d: %d, %.3f%%" $f $nc $qnc]
    puts $out  [format "Frame %d: %d, %.3f%%" $f $nc $qnc]
}
close $out
