#Half-assed diffraction calculation.

#Starters - calculate the location of the incoming wavefront.
#How to do this?  Trace a ray.  Compute slopes between surface 0 (located
#outside the telescope) and surface 1 (a plane at the front of the mirror).

#Note: I have 3 methods below that all do the same thing (and get the same
#answer, hopefully!

#################################################################
#I now need to run "ray" externally.
proc waveFront {optic xmm ymm xfract yfract icolor} {

#Trace ray
   ray $optic $xmm $ymm $xfract $yfract $icolor 0

#Get positions
   set xc [exprGet $optic.diagram->xray<0>]
   set yc [exprGet $optic.diagram->yray<0>]
   set zc [exprGet $optic.diagram->zray<0>]
   set xs [exprGet $optic.diagram->xray<1>]
   set ys [exprGet $optic.diagram->yray<1>]
   set zs [exprGet $optic.diagram->zray<1>]

#Compute slopes
   set mx [expr $xs - $xc]
   set my [expr $ys - $yc]
   set mz [expr $zs - $zc]
   set vnorm [expr sqrt($mx*$mx + $my*$my + $mz*$mz)]
   if {$vnorm == 0} return
   set mx [expr $mx/$vnorm]
   set my [expr $my/$vnorm]
   set mz [expr $mz/$vnorm]

#Compute point on wavefront that is at 0 phase.
#If zs is big negative, s is positive.
   set s [expr -($mx*$xs + $my*$ys + $mz*$zs)]

#echo s [format %.2f $s] mx [format %.2f $mx] my [format %.2f $my] \
	mz [format %.2f $mz]
   set xw [expr $xs + ($s-10000)*$mx]
   set yw [expr $ys + ($s-10000)*$my]
   set zw [expr $zs + ($s-10000)*$mz]

#Show I have setup xray<0>, etc, so that it is already at the same phase
#along the wavefront?  This would be easy to do in ray.
#echo xw = [format %.3f $xw] yw = [format %.3f $yw] zw = [format %.3f $zw]

#I will reset xray<0>, etc so it is done.
   handleSet $optic.diagram->xray<0> $xw
   handleSet $optic.diagram->yray<0> $yw
   handleSet $optic.diagram->zray<0> $zw
   return
   }

#########################################################
#Trace a bunch of rays and look for wavefront errors

proc waveErr {optic xmm color} {

#First, trace chief ray and get path length.
   set xsize [showFocal $optic $color xsize]
   
#Compute exit pupil
   opticInfo $optic $color

   set ymm 0
   waveFront $optic $xmm $ymm 0 0 $color
   set n [exprGet $optic.diagram->np]
   set n1 [expr $n-1]
   set n2 [expr $n-2]
   set xchief [exprGet $optic.diagram->xray<$n1>]
   set ychief [exprGet $optic.diagram->yray<$n1>]
   set zchief [exprGet $optic.diagram->zray<$n1>]
   set ep [showFocal $optic $color exit]

   set r2 [expr pow($zchief-$ep,2) + pow($xchief,2) + pow($ychief,2)]

#Compute path length
   set len 0
   set geomlen 0.
   set x2 [exprGet $optic.diagram->xray<0>]
   set y2 [exprGet $optic.diagram->yray<0>]
   set z2 [exprGet $optic.diagram->zray<0>]
   loop i 1 $n {
	set x1 $x2
	set y1 $y2
	set z1 $z2
	set x2 [exprGet $optic.diagram->xray<$i>]
	set y2 [exprGet $optic.diagram->yray<$i>]
	set z2 [exprGet $optic.diagram->zray<$i>]
	set surf1 [exprGet $optic.diagram->indx<[expr $i-1]>]
	set index [exprGet $optic.optic<$surf1>->n<$color>]
	set len [expr $len + sqrt(pow($x2-$x1,2) + pow($y2-$y1,2) + \
	   pow($z2-$z1,2))*[abs $index]]
	set geomlen [expr $geomlen + sqrt(pow($x2-$x1,2) + pow($y2-$y1,2) + \
	   pow($z2-$z1,2))]
	}

#index contains refraction index for last step of ray
   set lastindex $index
echo lastindex $lastindex
   set chieflen [expr $len - sqrt($r2*$lastindex*$lastindex)]

#2 directions
   plotInit a
   set xpoint ""
   set ypoint ""
   set scindex ""
   foreach dir "0. 1." {
      set sci [expr round($dir+2)]
      foreach fract "-1 -.8 -.6 -.4 -.2 0 .2 .4 .6 .8 1" {
	waveFront $optic $xmm $ymm [expr $fract*(1.-$dir)] \
	   [expr $fract*$dir] $color

#Position of off-chief ray in focal plane
	set xf [exprGet $optic.diagram->xray<$n1>]
	set yf [exprGet $optic.diagram->yray<$n1>]
	set zf [exprGet $optic.diagram->zray<$n1>]

#Slope at focal plane
	set mx [exprGet $optic.diagram->slope<0>]
	set my [exprGet $optic.diagram->slope<1>]
	set mz [exprGet $optic.diagram->slope<2>]

#Differences
	set dx [expr $xf -$xchief]
	set dy [expr $yf -$ychief]
	set dz [expr $zf -$zchief]
	set a 1.

#The sign in the following depends on whether the exit pupil position
#relative to the focal plane.  Better yet, use index of refraction
#to indicate direction.
	set sign [expr -1.*$lastindex]
###	if {$ep > $zchief} {set sign 1.} else {set sign -1.}
	set b [expr $sign*2.*($dx*$mx + $dy*$my + $dz*$mz)]
	set c [expr $dx*$dx + $dy*$dy + $dz*$dz - ($r2*$lastindex*$lastindex)]

#Distance from focal plane point to reference sphere along ray path.
	set s [expr (-$b + sqrt($b*$b - 4.*$a*$c))/(2.*$a)]

#Compute path length
	set len 0.
	set geomlen 0.
	set x2 [exprGet $optic.diagram->xray<0>]
	set y2 [exprGet $optic.diagram->yray<0>]
	set z2 [exprGet $optic.diagram->zray<0>]
	loop i 1 $n {
	   set x1 $x2
	   set y1 $y2
	   set z1 $z2
	   set x2 [exprGet $optic.diagram->xray<$i>]
	   set y2 [exprGet $optic.diagram->yray<$i>]
	   set z2 [exprGet $optic.diagram->zray<$i>]
	   set surf1 [exprGet $optic.diagram->indx<[expr $i-1]>]
	   set index [exprGet $optic.optic<$surf1>->n<$color>]
	   set len [expr $len + sqrt(pow($x2-$x1,2) + pow($y2-$y1,2) + \
		pow($z2-$z1,2))*[abs $index]]
	   set geomlen [expr $geomlen + sqrt(pow($x2-$x1,2) + pow($y2-$y1,2) \
		+ pow($z2-$z1,2))]
	   }

	set len [expr $len - $s]

	set err [expr $len - $chieflen]

#Error in waves
	set waveErr [format %10.4f \
	   [expr $err/([showWave $optic $color]*1.e-3)]]
	echo Dir $dir Fract [format %6.2f $fract], wavefront error [format \
	   %8.4f $err], in waves [format %8.2f $waveErr]
	lappend xpoint $fract
	lappend ypoint $waveErr
	lappend scindex $sci
	}
      }
   set vals [lsort -real $ypoint]
   set ymin [lindex $vals 0]
   set ymax [lindex $vals end]
   set max [expr 1.05*[max [abs $ymin] [abs $ymax]] + .05]
   pgSci 1
   pgEnv -1.05 1.1 -$max $max 0 0
   loop i 0 [llength $xpoint] {
	pgSci [lindex $scindex $i]
	pgPoint [lindex $xpoint $i] [lindex $ypoint $i] 3
	}
   return
   }

#########################################################
#Repeat above using Zernikes and zernikeWaveErr

proc waveErrZern {optic xmm color} {

set zern [genericNew ZERNIKE]
set thresh .01

#2 directions
   plotInit a
   set xpoint ""
   set ypoint ""
   set scindex ""
   set ymm 0
   set NORDER 2
   while {1} {
	zernikeFit $optic $zern $NORDER $xmm $ymm $color
	set err [exprGet $zern.fitErr]
	if {$err < $thresh} break
	if {$NORDER >= 14} break
	incr NORDER 2
	}

   foreach dir "0. 1." {
      set sci [expr round($dir+2)]
      foreach fract "-1 -.8 -.6 -.4 -.2 0 .2 .4 .6 .8 1" {
	set xfract [expr $fract*(1.-$dir)]
	set yfract [expr $fract*$dir]
#Error in waves
	set waveErr [format %10.4f [zernikeWaveErr $zern $xfract $yfract]]
	echo Dir $dir Fract [format %6.2f $fract], wavefront error in waves \
	   [format %8.2f $waveErr]
	lappend xpoint $fract
	lappend ypoint $waveErr
	lappend scindex $sci
	}
      }
   set vals [lsort -real $ypoint]
   set ymin [lindex $vals 0]
   set ymax [lindex $vals end]
   set max [expr 1.05*[max [abs $ymin] [abs $ymax]] + .05]
   pgSci 1
   pgEnv -1.05 1.1 -$max $max 0 0
   loop i 0 [llength $xpoint] {
	pgSci [lindex $scindex $i]
	pgPoint [lindex $xpoint $i] [lindex $ypoint $i] 3
	}
   zernikeZero $zern
   genericDel $zern
   return
   }

#########################################################
#Repeat above using C calculation of wave front errors.

proc waveErrChief {optic xmm color} {

set wave [genericNew WAVEFRONT]

#2 directions
   plotInit a
   set xpoint ""
   set ypoint ""
   set scindex ""
   set ymm 0
   chiefRay $optic $wave $xmm $ymm $color

   foreach dir "0. 1." {
      set sci [expr round($dir+2)]
      foreach fract "-1 -.8 -.6 -.4 -.2 0 .2 .4 .6 .8 1" {
	set xfract [expr $fract*(1.-$dir)]
	set yfract [expr $fract*$dir]
#Error in waves
	set waveErr [format %10.4f [offRay $optic $wave $xfract $yfract 0]]
	echo Dir $dir Fract [format %6.2f $fract], wavefront error in waves \
	   [format %8.2f $waveErr]
	lappend xpoint $fract
	lappend ypoint $waveErr
	lappend scindex $sci
	}
      }
   set vals [lsort -real $ypoint]
   set ymin [lindex $vals 0]
   set ymax [lindex $vals end]
   set max [expr 1.05*[max [abs $ymin] [abs $ymax]] + .05]
   pgSci 1
   pgEnv -1.05 1.1 -$max $max 0 0
   loop i 0 [llength $xpoint] {
	pgSci [lindex $scindex $i]
	pgPoint [lindex $xpoint $i] [lindex $ypoint $i] 3
	}
   genericDel $wave
   return
   }

########################################################################
#Encircled energy for CASS design. (actually, 3.5m )

proc ee {pix} {
   global PIXMM
   set fr 10.102
   set wave .5
   set x [expr $pix*$PIXMM]
   set v [expr 3.1416*$x/($wave*1.e-3*$fr)]
   set j0 [exec bes0 $v]
   set j1 [exec bes1 $v]
   set ee [expr 1. - $j0*$j0 - $j1*$j1]
   if {$v == 0} {
	set U 1.
   } else {
	set U [expr 2.*$j1/$v]
	}

#Central intensity
   set A [expr 1000.*3.1416*pow($PIXMM/(2.*$wave*1.e-3*$fr),2)]
   set psf [expr $A*pow($U,2)]
   echo x = [format %.3f $x] v/pi = [format %.3f [expr $v/3.1416]] \
	ee =  [format %.3f $ee] psf = [format %.4f $psf]
   }
