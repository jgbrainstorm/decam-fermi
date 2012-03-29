#######################################################################
#Differential refraction
#Compute diff. refr; input is 2 wavelengths in microns; zenith angle in degrees

proc diffRefract {w1 w2 z} {
   set d1 [expr [refract $w1] *tan($z/57.3)]
   set d2 [expr [refract $w2] *tan($z/57.3)]
   set diff [expr ($d2 - $d1) * 206265.]
   return [format %.2f $diff]
   }

#######################################################################
#n-1 for w in microns.  Corrected for ~7000 feet.  The correction factor is
#a guess.
proc refract {w} {
   return [expr (50./58.)*(2.735e-4 + 131.4182/pow($w*1.e4,2) + \
	2.76249e8/pow($w*1.e4,4))]
   }


