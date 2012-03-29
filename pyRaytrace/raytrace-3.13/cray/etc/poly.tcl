#Evaluate z height of a surface from curvature, a2, a4, a6, and a8.

proc poly {hndl surf rad} {
   set a2 [showparam $hndl $surf 8]
   set a4 [showparam $hndl $surf 9]
   set a6 [showparam $hndl $surf 10]
   set a8 [showparam $hndl $surf 11]
   set r2 [expr $a2*[pow $rad 2]]
   set r4 [expr $a4*[pow $rad 4]]
   set r6 [expr $a6*[pow $rad 6]]
   set r8 [expr $a8*[pow $rad 8]]
#Curvature
   set c [showparam $hndl $surf 1]
   if {$c != 0} {
	set R [expr 1./$c]
	if {$rad < $R} {
	   set z [expr $R - sqrt($R*$R - $rad*$rad)]
	   if {$c <0.} {set z [expr -1*$z]}
	} else {set z 0}
   } else {set z 0}
   return [expr $z + $r2 + $r4 + $r6 + $r8]
   }

#####################################################################
#Evaluate z height of a surface including all of above plus the real
#z height.
#This is used to help evaluate longitudinal chromatic aberration.

proc zpoly {hndl surf rad} {
   set zpos [showparam $hndl $surf 5]
   return [format %.3f [expr $zpos + [poly $hndl $surf $rad]]]
   }

##################################################################
#A real hack to get focal plane height for different filters
#Wire in the actual solutions.

proc zeval {hndl rad} {
   set list ""
   foreach ifil "2 3 4 5" {
     if {$ifil == 2} {
	setparam $hndl 10 5 756.49
	setparam $hndl 10 8 -4.73e-6
	setparam $hndl 10 9 -2.62e-10
	setparam $hndl 10 10 6.77e-16
     } elseif {$ifil == 3} {
	setparam $hndl 10 5 756.60
	setparam $hndl 10 8 -8.78e-6
	setparam $hndl 10 9 -2.54e-10
	setparam $hndl 10 10 6.35e-16
     } elseif {$ifil == 4} {
	setparam $hndl 10 5 756.67
	setparam $hndl 10 8 -1.14e-5
	setparam $hndl 10 9 -2.48e-10
	setparam $hndl 10 10 6.09e-16
     } elseif {$ifil == 5} {
	setparam $hndl 10 5 756.78
	setparam $hndl 10 8 -1.52e-5
	setparam $hndl 10 9 -2.40e-10
	setparam $hndl 10 10 5.71e-16
     } else {
	echo Unsupported filter $ifil
	return
	}
     set z [zpoly $hndl 10 $rad]
     lappend list $z
     }
   set min [eval min $list]
   set max [eval max $list]
   set diff [expr abs($max-$min)]
   return [format %0.3f $diff]
   }

