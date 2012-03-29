#Code to read baffle data (e.g., baffle.dat) and populate some data structures.

proc baffleRead {file} {
   set fid [open $file]
   set chain [chainNew BAFFLE]
   while {1} {
	set line [string trim [gets $fid]]
	if {[eof $fid]} break
	if {[string length $line] == 0} continue
	if {[string index $line 0] == "#"} continue
	set x1 [lindex $line 0]
	set y1 [lindex $line 1]
	set x2 [lindex $line 2]
	set y2 [lindex $line 3]
	set r1 [expr $x1*25.4]
	set r2 [expr $x2*25.4]
	set z1 [expr -1*$y1*25.4]
	set z2 [expr -1*$y2*25.4]
 	set hndl [genericNew BAFFLE]

#Is it an annulus?
	if {$z1 == $z2} {
	   handleSet $hndl.type ANNULUS
	   handleSet $hndl.r1 [min $r1 $r2]
	   handleSet $hndl.r2 [max $r1 $r2]
	   handleSet $hndl.z $z1
	} else {
	   set m [expr ($r2-$r1)/($z2-$z1)]
	   set z0 [expr -($z2*$r1-$z1*$r2)/($m*($z2-$z1))]
	   handleSet $hndl.type CONE
	   handleSet $hndl.m $m
	   handleSet $hndl.z0 $z0
	   handleSet $hndl.z1 [min $z1 $z2]
	   handleSet $hndl.z2 [max $z1 $z2]
	   }
	chainAdd $chain $hndl
	}
   close $fid
   return $chain
   }
