# q-value separato per classe di separazione in sequenza.
# Da usare come nativecontacts.tcl (molecola gia' caricata, es. via source da analysisMD.tcl).
# Scrive q_split.txt:  frame  n_tot q_tot  n_loc q_loc  n_med q_med  n_ter q_ter
#   loc = separazione 3-4   (struttura secondaria locale)
#   med = separazione 5-9
#   ter = separazione >=10  (impacchettamento terziario: il segnale piu' sensibile)

set cutoff 7.0

set protein [atomselect top "protein and name CA"]
foreach i [$protein get index] r [$protein get resid] { set resid($i) $r }

proc classe {sep} {
    if {$sep < 3}   { return "" }
    if {$sep <= 4}  { return loc }
    if {$sep <= 9}  { return med }
    return ter
}

proc contactSets {sel cutoff} {
    global resid
    set c [measure contacts $cutoff $sel]
    set out [dict create loc {} med {} ter {}]
    foreach i [lindex $c 0] j [lindex $c 1] {
        set k [classe [expr {abs($resid($i) - $resid($j))}]]
        if {$k eq ""} continue
        set p [expr {$i < $j ? "$i $j" : "$j $i"}]
        dict lappend out $k $p
    }
    foreach k {loc med ter} { dict set out $k [lsort -unique [dict get $out $k]] }
    return $out
}

set nf [molinfo top get numframes]

$protein frame 0
set ref [contactSets $protein $cutoff]
foreach k {loc med ter} {
    set n($k) [llength [dict get $ref $k]]
    foreach p [dict get $ref $k] { set nativo($k,$p) 1 }
}
set n(tot) [expr {$n(loc) + $n(med) + $n(ter)}]

set out [open ./q_split.txt w]
puts       "Native: tot=$n(tot)  loc=$n(loc)  med=$n(med)  ter=$n(ter)"
puts $out "# Native: tot=$n(tot)  loc=$n(loc)  med=$n(med)  ter=$n(ter)"
puts $out "# frame n_tot q_tot n_loc q_loc n_med q_med n_ter q_ter"

for {set f 0} {$f < $nf} {incr f} {
    $protein frame $f
    set cur [contactSets $protein $cutoff]
    set tot 0
    foreach k {loc med ter} {
        set c($k) 0
        foreach p [dict get $cur $k] {
            if {[info exists nativo($k,$p)]} { incr c($k) }
        }
        incr tot $c($k)
    }
    set riga [format "%d %d %.3f" $f $tot [expr {100.0*$tot/$n(tot)}]]
    foreach k {loc med ter} {
        append riga [format " %d %.3f" $c($k) [expr {$n($k) > 0 ? 100.0*$c($k)/$n($k) : 0}]]
    }
    puts $out $riga
    if {$f % 100 == 0} { puts "frame $f: $riga" }
}
close $out
