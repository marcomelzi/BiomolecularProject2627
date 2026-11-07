# Version Jan-2022

# ---------------- PARAMETRI DA CONTROLLARE ----------------
set psf   ionized.psf    ;# deve essere lo STESSO file della riga "structure" nel tuo .conf
set dcd   prot3.dcd      ;# traiettoria da analizzare
set step  1              ;# analizza 1 frame ogni $step (metti 10 se la dcd e' grossa)
# ----------------------------------------------------------

mol new $psf type psf waitfor all
mol addfile $dcd type dcd waitfor all

# numero di frame effettivamente caricati (MANCAVA nello script originale)
set nf [molinfo top get numframes]
puts "Frame caricati: $nf"
if {$nf == 0} {
    puts "ERRORE: nessun frame letto da $dcd  (psf/dcd sbagliati o percorso errato)"
    exit
}

set all [atomselect top all]
set bb [atomselect top backbone]
set bb0 [atomselect top backbone]
set protein [atomselect top protein]
set ca [atomselect top "name CA"]
set numres [$ca num]


$all frame 0

# Align trajectory to first frame
puts "Align"
set reference [atomselect top "backbone" frame 0]
set compare [atomselect top "backbone"]
for {set i 0} {$i < $nf} {incr i} {
        $all frame $i
        $compare frame $i
        # compute the transformation    
        set trans_mat [measure fit $compare $reference]
        $all move $trans_mat
}
  
#### Calculate RMSD
if {1} {
set out [open rmsd.dat w]
puts "Calculate RMSD"
	for {set i 0} {$i < $nf} {incr i $step} {
		$bb frame $i	
		$bb0 frame 0	
		set rmsd [measure rmsd $bb0 $bb]
		#puts $i	
		puts $out "$i $rmsd"
	}
close $out
}


#### Calculate RMSF
if {1} {
set out [open rmsfCA.dat w]
puts "Calculate RMSF"
set rmsfCA [measure rmsf $ca]
set firstAA [lindex [$ca get resid] 0]
for {set i 0} {$i < [$ca num]} {incr i} {
	puts $out "[expr $i+$firstAA] [lindex $rmsfCA $i]"
	}
close $out
}


#### Calculate Hbonds
if {0} {
package require hbonds
puts "Calculate Hbonds"
hbonds -sel1 $bb -upsel no -writefile yes -polar yes -plot no -frames 0:$step:[expr $nf-1]
}


#### Calculate SASA
if {0} {
set out [open sasa.dat w]
puts "Calculate SASA"
	for {set i 0} {$i < $nf} {incr i $step} {
		$protein frame $i	
		set sasa [measure sasa 1.4 $protein]
		#puts $i	
		puts $out "$i $sasa"
	}
close $out
}


#### Calculate Rgyr
if {0} {
set out [open rgyr.dat w]
puts "Rgyr"
for {set i 0} {$i < $nf} {incr i $step} {
	$protein frame $i
	set rgyr [measure rgyr $protein]
	#puts $i	
	puts $out "$i $rgyr"
}
close $out
}

#### Calculate Beta sheet
if {0} {
set out [open beta.dat w]
puts "Calculate Beta"
for {set i 0} {$i < $nf} {incr i $step} {
        animate goto $i
        mol ssrecalc top
        set st [$ca get structure]
        set Beta [regexp -all -- E $st]
        #Beta: % of beta over residues
        set BetaPerc [ expr 1.00 * $Beta / $numres ]
        #puts $i 
        puts $out "$i $BetaPerc $Beta"
}
close $out
}


#### Calculate Alpha
if {0} {
set out [open alpha.dat w]
puts "Calculate Alpha"
for {set i 0} {$i < $nf} {incr i $step} {
        animate goto $i
        mol ssrecalc top
        set st [$ca get structure]
        set Alpha [regexp -all -- H $st]
        #Alphaperc: % of alpha over residues
        set AlphaPerc [ expr 1.00 * $Alpha / $numres ]
        #puts $i 
        puts $out "$i $AlphaPerc $Alpha"
}
close $out
}

exit




