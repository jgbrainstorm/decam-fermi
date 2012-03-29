#Procedure to check for vignetting due to baffles

proc rayCheck {optic baffle pos filter} {
   set finner [exprGet $optic.tel->finner]
   set npt 0
   set ntot 0
   loop i 0 101 {
	set fract [format %.2f [expr $i/100.]]
	echo $fract
	if {abs($fract) < $finner} continue
	set nj [expr $fract*314]
	loop j 0 $nj {
	   set ang [expr $j*3.14/$nj]
	   set xfract [expr $fract*cos($ang)]
	   set yfract [expr $fract*sin($ang)]
	   ray $optic $pos 0 $xfract $yfract $filter
	   set line [baffleCheck $optic $baffle]
	   incr ntot
	   if {$line == ""} {
		incr npt
		}
	   }
	}
   return [expr double($npt)/$ntot]
   }
