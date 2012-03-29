#NOTE: This file is NOT loaded at startup.  These are just TEST routines.
#
# Some tcl procs for experimenting with zernikes

# Wavefront error is sum of terms.  We need derivatives in r, theta in order
# to compute gradient in Cartesian coords.

proc dWdRdT {Zin r theta} {
   upvar 1 $Zin Z
   set cs [cos $theta]
   set sn [sin $theta]
   set cs2 [cos 2.*$theta]
   set sn2 [sin 2.*$theta]
   set r2 [expr $r*$r]
   set dWdR [expr $Z(1)*$cs + $Z(2)*$sn + $Z(3)*4.*$r + $Z(4)*2.*$r*$cs2 \
	+ $Z(5)*2.*$r*$sn2 + $Z(6)*(9.*$r2-2.)*$cs + $Z(7)*(9.*$r2-2.)*$sn]
   set dWdT [expr -$Z(1)*$sn + $Z(2)*$cs - $Z(4)*$r2*2.*$sn2 + \
	$Z(5)*$r2*2.*$cs2 - $Z(6)*(3.*$r2-2.)*$r*$sn + \
	$Z(7)*(3.*$r2-2.)*$r*$cs]
   return "$dWdR $dWdT"
   }

#Compute gradient vector

proc zGrad {Zin r theta} {
   upvar 1 $Zin Z
   set x [expr $r*cos($theta)]
   set y [expr $r*sin($theta)]
   set r2 [expr $r*$r]
   set grads [dWdRdT Z $r $theta]
   set dWdR [lindex $grads 0]
   set dWdT [lindex $grads 1]
   if {$r == 0} {return "0 0"}
   set xGrad [expr $dWdR*$x/$r - $dWdT*$y/$r2]
   set yGrad [expr $dWdR*$y/$r + $dWdT*$x/$r2]
   return "$xGrad $yGrad"
   }

#Top level proc to plot up a spot diagram for a zernike polynimial
#Must specify Z(0) thru Z(7) in a global variable

proc zernike {} {
   global Z
   pgEnv -2 2 -2 2 1 0
   loop i 3 10 {
	set r [expr ($i+1.)/10.]
	loop j 0 30 {
	   set ang [expr $j*(2.*3.141593)/30.]
	   set grads [zGrad Z $r $ang]
	   set x [expr $r*cos($ang) - [lindex $grads 0]]
	   set y [expr $r*sin($ang) - [lindex $grads 1]]
	   pgPoint $x $y 1
	   }
	}
   return
   }

######################################################################
#List m, n, N, j values that are used in zernike.c

proc zlist {NORDER} {
   set M [expr $NORDER/2]
   for {set m 0} {$m <= $M} {incr m} {
	for {set n $m} {$n <= $NORDER-$m} {set n [expr $n+2]} {
	   set k [expr ($n-$m)/2]
	   set N [expr ($n+$m)]
	   set j [expr ($N/2)*($N/2) + 2*$k]
#	   echo N $N m $m n $n j $j
	   set line($j)  "N $N m $m n $n"
	   if {$m > 0} {
		incr j
#		echo N $N m $m n $n j $j
		set line($j) "N $N m $m n $n"
		}
	   }
	}
   foreach j [lsort -integer [array names line]] {
	echo j $j $line($j)
	}
   return
   }

