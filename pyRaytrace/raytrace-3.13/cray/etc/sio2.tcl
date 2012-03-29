proc sio2 {n} {
   set xlist {.3650   .4047   .4678   .5086   .5770   .6563   .7065   .8521   \
	.8943}
   set ylist {1.47454 1.46962 1.46429 1.46182 1.46008 1.45637 1.45514 1.45246 \
	1.45183}
   set coeff [polyfit $xlist $ylist $n]
   set res [polycomp $xlist $ylist $coeff]
   set rms [polyrms $res]
   echo rms is $rms
   set waves [list 1.01398 0.85211 0.70652 0.65627 0.58756 \
     0.54607 0.48613 0.43583 0.40466 0.36501]
   set index ""
   foreach wave $waves {
	set val [polyeval $wave $coeff]
	lappend index $val
	}
   echo [eval format %6.4f%7.4f%7.4f%7.4f%7.4f%7.4f%7.4f%7.4f%7.4f%7.4f \
	$waves]
   echo [eval format %6.4f%7.4f%7.4f%7.4f%7.4f%7.4f%7.4f%7.4f%7.4f%7.4f \
	$index]
   }
