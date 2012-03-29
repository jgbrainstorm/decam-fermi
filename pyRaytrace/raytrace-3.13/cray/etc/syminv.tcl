#TCL version of least squares fitting.

########################################################################
#This shows how to call syminv
proc syminvtest {} {

#Number of parameters
   set nparam 3
   set l 0
   loop j 0 $nparam {
	set rhs($j) 0.
	for {set k 0} {$k <= $j} {incr k} {
	   set cmat($l) 0.
	   incr l
	   }
	}

#Coeffs of polynomials
   set a 1.
   set b 3.
   set c -2.

#Loop through observables.
   set nstep 2
   for {set i -$nstep} {$i <= $nstep} {incr i} {
	set x [expr ($i*1.)/$nstep]

#Zero out derivatives
	loop j 0 $nparam {
	   set deriv($j) 0.
	   }

#O-C and derivs.
	set obs [expr $a + $b*$x + $c*pow($x,2)]
	set comp 0.
	set res [expr $obs-$comp]
	set l 0
	loop j 0 $nparam {
	   set deriv($j) [expr pow($x,$j)]
	   set rhs($j) [expr $rhs($j) + $res * $deriv($j)]
	   for {set k 0} {$k <= $j} {incr k} {
		set cmat($l) [expr $cmat($l) + $deriv($j) * $deriv($k)]
		incr l
		}
	   }
	}
parray cmat
parray rhs
   syminv cmat rhs
   loop j 0 $nparam {
	echo j $j rhs $rhs($j)
	}
   return
   }
########################################################################
#Invert matrix

proc syminv {_incmat _inrhs} {
   upvar 1 $_incmat cmat
   upvar 1 $_inrhs rhs
   set dim [array size rhs]
echo dim is $dim
   loop i 0 $dim {
	set use($i) 0
	}

   set pvrwb 0.
   set pivot 0.
   loop n 0 $dim {
	set dmax 0.
	set j 0
	for {set i 0} {$i < $dim} {incr i; incr j [expr $i+1]} {
	   set diag [expr abs($cmat($j))]
	   if {$diag-$dmax <= 0.} continue
	   if {$use($i) > 0} continue
	   set dmax $diag
	   set k $i
	   set pivot $cmat($j)
	   }
	if {$dmax <= 0.} {return 0}
	set use($k) 1
	set l [expr $k*($k+1)/2]
	loop i 0 $dim {
	   if {$i-$k < 0} {
		set pvrow($i) $cmat($l)
		set pvcol($i) [expr $cmat($l)/$pivot]
		set cmat($l) 0.
		incr l
		if {$use($i) == 0} continue
		set pvcol($i) [expr -1.*$pvcol($i)]
		continue
		}
	   if {$i-$k == 0} {
		set pvrow($i) 1.
		set pvcol($i) [expr -1./$pivot]
		set cmat($l) 0.
		set pvrwb $rhs($i)
		incr l [expr $i+1]
		set rhs($i) 0.
		continue
		}
	   set pvrow($i) $cmat($l)
	   set pvcol($i) [expr $cmat($l)/$pivot]
	   set cmat($l) 0.
	   if {$use($i) > 0} {
		set pvrow($i) [expr -1.*$pvrow($i)]
		}
	   incr l [expr $i+1]
	   }
	set m 0
	loop i 0 $dim {
	   set rhs($i) [expr $rhs($i) - $pvrwb*$pvcol($i)]
	   for {set h 0} {$h <= $i} {incr h} {
		set cmat($m) [expr $cmat($m) - $pvcol($i)*$pvrow($h)]
		incr m
		}
	  }
	}
   return 1
   }


