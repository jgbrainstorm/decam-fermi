#Compute dispersion for ZnSe, SIO2 and design a prism that has constant
#R(theta) - i.e., delta-lambda/lambda is constant over 1-2 microns.

#Compute dispersion function
proc jdemDisp {} {
   global sio2 znse
   if {[info exists sio2]} {unset sio2}
   if {[info exists znse]} {unset znse}
   set sio2() ""
   set znse() ""


   set xlist ""
   set ylist ""
   loop i 0 11 {
	set wave [expr 1.+($i/10.)]
	set dwave [expr $wave+.02]
	lappend sio2() $wave
	lappend znse() $wave
	set sio2($wave) [expr ([glass SIO2 $dwave]-[glass SIO2 $wave])/.02]
	set znse($wave) [expr ([glass ZnSe $dwave]-[glass ZnSe $wave])/.02]
	}
   return
   }

#####################################################################
proc prismFit {} {
   global sio2 znse

#Model is a*sio2 + b*znse + c = 1/lambda
   set l 0
   loop j 0 3 {
	set rhs($j) 0.
	for {set k 0} {$k <= $j} {incr k} {
	   set cmat($l) 0
	   incr l
	   }
	}
   foreach wave $sio2() {
	set obs($wave) [expr -1./$wave]
	set prt(0) $sio2($wave)
	set prt(1) $znse($wave)
	set prt(2) 1.
	set l 0
	loop j 0 3 {
	   set rhs($j) [expr $rhs($j) + $prt($j)*$obs($wave)]
	   for {set k 0} {$k <= $j} {incr k} {
		set cmat($l) [expr $cmat($l) + $prt($j)*$prt($k)]
		incr l
		}
	   }
	}
   syminv cmat rhs
   echo a = $rhs(0) b = $rhs(1)

   foreach wave $sio2() {
	set comp [expr $rhs(0)*$sio2($wave) + $rhs(1)*$znse($wave) + $rhs(2)]
	set res [expr $obs($wave) - $comp]
	echo [format %.1f $wave] obs [format %.4f $obs($wave)] \
	    comp [format %.1f $comp]
	}

   return
   }
