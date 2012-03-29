#Compute angular separation given ra, dec in degrees

proc plateDiff {ra1 dec1 ra2 dec2} {
   set pi 3.1415926536

   set r1 [expr ($ra1)*$pi/180.]
   set d1 [expr ($dec1)*$pi/180.]
   set r2 [expr ($ra2)*$pi/180.]
   set d2 [expr ($dec2)*$pi/180.]

   set x1 [expr cos($r1)*cos($d1)]
   set y1 [expr sin($r1)*cos($d1)]
   set z1 [expr sin($d1)]
   set x2 [expr cos($r2)*cos($d2)]
   set y2 [expr sin($r2)*cos($d2)]
   set z2 [expr sin($d2)]
   set xp [expr $y1*$z2 - $z1*$y2]
   set yp [expr $z1*$x2 - $x1*$z2]
   set zp [expr $x1*$y2 - $x2*$y1]
   set sep [expr sqrt(pow($xp,2) + pow($yp,2) + pow($zp,2))]
   set sep [expr $sep*180./$pi]
   return $sep
   }
