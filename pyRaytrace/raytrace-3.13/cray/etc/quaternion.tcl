#Procedures to manipulate quaternions and Euler angles

#euler2quat	omega incl phi
#euler2quatdeg	omegad incld phid
#quat2euler	q
#quat2eulerdeg	q
#quatprod	q1 q2
#quatconj	q
#quat2rot	q
#rotat		vec c

#Convert Euler angles to a quaternion
#Omega is longitude of ascending node.
#incl is inclination
#phi is argument of periapse or equivalent

#If I do Euler angles in succession, then
#R(total) = R(Omega) * R(incl) * R(phi)

#I can compute R(total) by calling quatprod twice:
#quatprod [quatprod <q1> <q2>] <q3>
#where q1 = euler2quat <Omega> 0. 0.  etc.

#Input in degrees or radians?  I guess I will use radians, since that is
#how angles are stored in the OPTIC structure.  But I will also provide
#a deg version as well.

proc euler2quat {omega incl phi} {
   set co [expr cos($omega/2.)]
   set so [expr sin($omega/2.)]
   set ci [expr cos($incl/2.)]
   set si [expr sin($incl/2.)]
   set cp [expr cos($phi/2.)]
   set sp [expr sin($phi/2.)]

   set a [expr  $co * $ci * $cp - $so * $ci * $sp]
   set b [expr  $co * $si * $cp + $so * $si * $sp]
   set c [expr -$co * $si * $sp + $so * $si * $cp]
   set d [expr  $co * $ci * $sp + $so * $ci * $cp]

   return [list $a $b $c $d]
   }

################################################

#A degrees version
proc euler2quatdeg {omegad incld phid} {
   set DEG2RAD [expr 3.14159265358979/180.]
   set omega [expr $omegad*$DEG2RAD]
   set incl [expr $incld*$DEG2RAD]
   set phi [expr $phid*$DEG2RAD]
   return [euler2quat $omega $incl $phi]
   }

################################################

#Input is a list
proc quat2euler {q} {
   set a [lindex $q 0]
   set b [lindex $q 1]
   set c [lindex $q 2]
   set d [lindex $q 3]

   set omega [format %.6f [expr atan2($a*$c + $b*$d, -$c*$d + $a*$b)]]
   set phi [format %.6f [expr atan2($b*$d - $a*$c, $a*$b + $c*$d)]]
   set ci [expr 1. - 2.*($b*$b + $c*$c)]
   if {$ci > 1.} {set ci 1.}
   if {$ci < -1.} {set ci -1.}
   set incl [format %.6f [expr acos($ci)]]
   return [list $omega $incl $phi]
   }

################################################

#Degree version of quat2euler
proc quat2eulerdeg {q} {
   set RAD2DEG [expr 180./3.14159265358979]
   set list [quat2euler $q]
   set omega [lindex $list 0]
   set incl [lindex $list 1]
   set phi [lindex $list 2]
   set omegad [format %.5f [expr $omega*$RAD2DEG]]
   set incld [format %.5f [expr $incl*$RAD2DEG]]
   set phid [format %.5f [expr $phi*$RAD2DEG]]
   return [list $omegad $incld $phid]
   }

################################################

#Multiply two quaternions
proc quatprod {q1 q2} {
   set a1 [lindex $q1 0]
   set b1 [lindex $q1 1]
   set c1 [lindex $q1 2]
   set d1 [lindex $q1 3]

   set a2 [lindex $q2 0]
   set b2 [lindex $q2 1]
   set c2 [lindex $q2 2]
   set d2 [lindex $q2 3]

   set a [expr $a1*$a2 - $b1*$b2 - $c1*$c2 - $d1*$d2]
   set b [expr $a1*$b2 + $b1*$a2 + $c1*$d2 - $d1*$c2]
   set c [expr $a1*$c2 - $b1*$d2 + $c1*$a2 + $d1*$b2]
   set d [expr $a1*$d2 + $b1*$c2 - $c1*$b2 + $d1*$a2]

   return [list $a $b $c $d]
   }

################################################

#Conjugate of a quaternion
proc quatconj {q} {
   set a [lindex $q 0]
   set b [lindex $q 1]
   set c [lindex $q 2]
   set d [lindex $q 3]

   set b [expr -1.*$b]
   set c [expr -1.*$c]
   set d [expr -1.*$d]

   return [list $a $b $c $d]
   }

###############################################

#Compute rotation matrix from a quaternion
proc quat2rot {q} {

   set a [lindex $q 0]
   set b [lindex $q 1]
   set c [lindex $q 2]
   set d [lindex $q 3]

#Matrix is a list of sublists.  Each sublist is a row in the matrix
#c12 means row 1, col 2 of matrix.
   set c11 [expr 1. - 2.*($c*$c + $d*$d)]
   set c12 [expr 2.*($b*$c - $a*$d)]
   set c13 [expr 2.*($a*$c + $b*$d)]
   set list1 [list $c11 $c12 $c13]

   set c21 [expr 2.*($b*$c + $a*$d)]
   set c22 [expr 1. - 2.*($b*$b + $d*$d)]
   set c23 [expr 2.*($c*$d - $a*$b)]
   set list2 [list $c21 $c22 $c23]

   set c31 [expr 2.*($b*$d - $a*$c)]
   set c32 [expr 2.*($a*$b + $c*$d)]
   set c33 [expr 1. - 2.*($b*$b + $c*$c)]
   set list3 [list $c31 $c32 $c33]

   return [list $list1 $list2 $list3]
   }

#####################################

#Apply rotation to a 3-vector.  Input is vector, rotation matrix
proc rotat {vec c} {
   set x [lindex $vec 0]
   set y [lindex $vec 1]
   set z [lindex $vec 2]

   set list1 [lindex $c 0]
   set list2 [lindex $c 1]
   set list3 [lindex $c 2]

   set c11 [lindex $list1 0]
   set c12 [lindex $list1 1]
   set c13 [lindex $list1 2]

   set c21 [lindex $list2 0]
   set c22 [lindex $list2 1]
   set c23 [lindex $list2 2]

   set c31 [lindex $list3 0]
   set c32 [lindex $list3 1]
   set c33 [lindex $list3 2]

   set xout [expr $c11*$x + $c12*$y + $c13*$z]
   set yout [expr $c21*$x + $c22*$y + $c23*$z]
   set zout [expr $c31*$x + $c32*$y + $c33*$z]

   return [list $xout $yout $zout]
   }
