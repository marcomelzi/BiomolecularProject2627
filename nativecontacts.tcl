## The atom selection
 set protein [ atomselect top "protein and name CA"]
 
 ## Use first trajectory frame as a reference
 $protein frame 0
 set ref [ prepareNativeContacts 7 $protein]
 
## Output file
set out [open ./q_value/contacts.txt w] 

 ## Get the number of native contacts (for computing their fraction)
 set nnc [ llength $ref]
 puts "There are $nnc contacts in the native state"
 puts $out "There are $nnc contacts in the native state"

               # Now, for each frame,
 forFrames fno $protein {
               # compute number of native contacts,
         set nc [measureNativeContacts $ref 7 $protein]  
               # their fraction,
         set qnc [ expr 100.0 * $nc / $nnc ]
               # and print both.
         puts [ format "Frame %d: %f, %.3f%%" $fno $nc $qnc ]
	 puts $out [ format "Frame %d: %f, %.3f%%" $fno $nc $qnc ]
 }

close $out

