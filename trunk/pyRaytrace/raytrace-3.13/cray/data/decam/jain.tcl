#More weak lensing code in support of Jarvis & Jain
#
#wplot <optic> ifil	Compute rms residuals for static pattern
#wplot <optic1> <optic2> ifil	Compute rms residuals for difference of 2
#				designs.
#wglobal <optic>	Compute static rms over entire camera
#wglobal <optic1> <optic2>	Copmute differential rms over entire camera.

#Compute whiskers on a grid.

proc whiskGrid {hndl ifil} {
   set xrad [showFocal $hndl $ifil xrad]
   set yrad [showFocal $hndl $ifil yrad]

#Limits of inscribed rectangle
   set xrect [expr $xrad/sqrt(2.)]
   set yrect [expr $yrad/sqrt(2.)]

   set fid [open whisker w]
   loop i 0 1024 {
	set outlist ""
	loop j 0 1024 {
	   set x [expr ($i-511.5)/512.*$xrect]
	   set y [expr ($j-511.5)/512.*$yrect]
	   set list [whisker $hndl $x $y $ifil]
	   lappend outlist [lindex $list 0] [lindex $list 1]
	   }
	puts $fid $outlist
	}
   close $fid
   return
   }

#######################################################################
#Divide focal plane into a mosaic of detectors.
#Design originally has 8 filters.  I will add 2 dummys, then start
#mosaic numbering from 11.

proc decamDesign {optic} {

#Reference filter
   set rfil 5

   set wave [showFocal $optic $rfil wave]
   set scale [showFocal $optic $rfil scale]
   set z [showSurf $optic 15 z]
   waveAdd $optic $wave $rfil
   waveAdd $optic $wave $rfil

#Add the individual detectors in filters 11 to 72
#Number of surfaces now - must fill in optical index for new filters

#   set focalid 15
   set surflist [surfGet $optic $rfil]
   set focalid [lindex $surflist end]
   set z [showSurf $optic $focalid z]

#Optical CCDs
   set ccdwidth [expr .015*4096.]
   set ccdheight [expr .015*2048.]

#Gaps between CCDs - this is approx 114 pixels.
   set xgap 2.45
   set ygap 3.096

#Include wavefront and guiding sensors - correct positions later
   set n(0) 7
   set n(1) 6
   set n(2) 6
   set n(3) 5
   set n(4) 6
   set n(5) 5
   set n(6) 2

#Row, col of 2Kx2K CCDs
   set 2k [list 4,0 4,5 5,0 5,4 6,0 6,1]

#2 sides of array
   set nfil 10
   set surfid $focalid
   foreach side "-1. 1." {

#Loop through 7 rows of CCDs
	loop j 0 7 {
	   set yc [expr (($ccdheight+$ygap)*($j+0.5))*$side]
	   loop i 0 $n($j) {
		set xc [expr (-($n($j)-1.)/2. + $i)*($ccdwidth + $xgap)]

#Is this a 2Kx2K CCD?  Adjust xc if so
		set width $ccdwidth
		set height $ccdheight
		if {[lsearch $2k $j,$i] >= 0} {
		   if {$xc < 0} {
			set xc [expr $xc + $ccdwidth/4.]
		  } else {
			set xc [expr $xc - $ccdwidth/4.]
			}
		   set width [expr $width/2.]
		   }
		incr nfil

#Use opticInsert to add a new surface.  Get new surfid.
		set newid [opticInsert $optic $surfid]
		setName $optic $newid FOCAL
		setSurf $optic $newid x $xc
		setSurf $optic $newid y $yc
		setSurf $optic $newid z $z
		waveAdd $optic $wave $rfil
		setIndex $optic $focalid $nfil 0
		setIndex $optic $newid $nfil [showIndex $optic $focalid $rfil]
		set surfid $newid

#Regression test - ncolor should be the same as nfil
		set ncolor [exprGet $optic.ncolor]
		if {$nfil != $ncolor} {
		   error "Gak! filter count off nfil = $nfil, ncolor = $ncolor"
		   }

#Position on sky.  Don't know if my signs are correct yet!
		set xoff [expr $xc*$scale/60.]
		set yoff [expr $yc*$scale/60.]
		setFocal $optic $nfil xoff $xoff
		setFocal $optic $nfil yoff $yoff
		setFocal $optic $nfil xsize [expr -$width/2.]
		setFocal $optic $nfil ysize [expr -$height/2.]
		setFocal $optic $nfil scale $scale
		setFocal $optic $nfil dist 0
		setFocal $optic $nfil rot 0
		opticInfo $optic $nfil

#Refine xoff, yoff
		ray $optic 0 0 0 0 $nfil 0
		set np [exprGet $optic.diagram->np]
		set np1 [expr $np-1]
		set xfocal [exprGet $optic.diagram->xray<$np1>]
		set yfocal [exprGet $optic.diagram->yray<$np1>]
		set dx [expr $xc-$xfocal]
		set dy [expr $yc-$yfocal]
		set dxoff [expr $dx*$scale/60.]
		set dyoff [expr $dy*$scale/60.]
		set xoff [expr $xoff + $dxoff]
		set yoff [expr $yoff + $dyoff]
		setFocal $optic $nfil xoff $xoff
		setFocal $optic $nfil yoff $yoff
		}
	   }
	}

   stopcomp $optic
   opticinc $optic 1
   return
   }

#########################################################################
#Plot <xx>-<yy> and <xy> vs. x position for any focal plane CCD.

proc wplot {optic args} {
   global SLEEP
   if {![info exists SLEEP]} {
	set SLEEP 0
	}
   echo Using global variable SLEEP timer of $SLEEP seconds.

   if {[llength $args] == 1} {
	set ifil [lindex $args 0]
	set noptic 1
   } elseif {[llength $args] == 2} {
	set optic2 [lindex $args 0]
	set ifil [lindex $args 1]
	set noptic 2
   } else {
	error "Usage: wplot <optic> \[<optic2>\] ifil"
	}
   global xx xy dil
   set xrad [showFocal $optic $ifil xrad]
   set yrad [showFocal $optic $ifil yrad]
   set n 20
   set ylist ""
   pgSci 1
   pgEnv [expr -1.*abs($xrad)] [expr abs($xrad)] -.15 .15 0 0
   for {set j -2} {$j <= 2} {incr j} {
	set y [expr abs($yrad)*$j/2.]
	lappend ylist $y
	set xlist ""
	set xxlist ""
	set xylist ""
	set dillist ""
	for {set i -$n} {$i <= $n} {incr i} {
	   set x [expr abs($xrad)*($i*1.)/$n]
	   lappend xlist $x
	   set list [polarize $optic $x $y $ifil]
	   set xx($x,$y) [lindex $list 0]
	   set xy($x,$y) [lindex $list 1]
	   set dil($x,$y) [lindex $list 2]

	   if {$noptic == 2} {
		set list [polarize $optic2 $x $y $ifil]
		set xx($x,$y) [expr $xx($x,$y) - [lindex $list 0]]
		set xy($x,$y) [expr $xy($x,$y) - [lindex $list 1]]
		set dil($x,$y) [expr $dil($x,$y) - [lindex $list 2]]
		}

	   lappend xxlist $xx($x,$y)
	   lappend xylist $xy($x,$y)

#e3 is actually square of dilution.  It works out that dil itself is more
#linear.  Oops, this does not work when differencing designs
	   if {$noptic == 1} {set dil($x,$y) [expr sqrt($dil($x,$y))]}
	   lappend dillist [expr $dil($x,$y)*.2]
	   }
	pgLine $xlist $xxlist
	pgSci 3
	pgLine $xlist $xylist
	pgSci 5
	pgLine $xlist $dillist
	pgSci 1
	}
   sleep $SLEEP
#Cache xlist and ylist info so I can factor out the bilinear fits to a
#separate proc
   set xx(xlist) $xlist
   set xx(ylist) $ylist
   set xy(xlist) $xlist
   set xy(ylist) $ylist
   set dil(noptic) $noptic
   set dil(xlist) $xlist
   set dil(ylist) $ylist
   set rms [wfit]
   return $rms
   }

############################################################

#I would like to make bilinear fit.
#I can collapse to 1d on each axis and do it that way.
proc wfit {} {
   global SLEEP
   if {![info exists SLEEP]} {
	set SLEEP 0
	}
   global xx xy dil

#Assume xlist and ylist are the same for all 3 arrays, which they are.
   set xlist $xx(xlist)
   set ylist $xx(ylist)
   set noptic $dil(noptic)
   set xrad [eval max $xlist]
   set yrad [eval max $ylist]

#Remove average from x, y to make them 0-centered
   set xavg 0.
   set yavg 0.
   set nx 0
   set ny 0
   foreach x $xlist {
	set xavg [expr $xavg + $x]
	incr nx
	}
   set xavg [expr $xavg/(1.*$nx)]

   foreach y $ylist {
	set yavg [expr $yavg + $y]
	incr ny
	}
   set yavg [expr $yavg/(1.*$ny)]

#Fit xx in x direction
   set rlist ""
   set obslist ""
   foreach x $xlist {
	lappend rlist [expr $x-$xavg]
	set sum 0.
	foreach y $ylist {
	   set sum [expr $sum + $xx($x,$y)]
	   }
	lappend obslist [expr $sum/(1.*$ny)]
	}
   set xcoeff [polyfit $rlist $obslist 2]
   set reslist [polycomp $rlist $obslist $xcoeff]
   set rms [polyrms $reslist]
   echo xx xrms $rms

#Cache xx coeff
   set xxxslope [lindex $xcoeff 1]

#Fit xx in y direction
   set rlist ""
   set obslist ""
   foreach y $ylist {
	lappend rlist [expr $y-$yavg]
	set sum 0.
	foreach x $xlist {
	   set sum [expr $sum + $xx($x,$y)]
	   }
	lappend obslist [expr $sum/(1.*$nx)]
	}
   set ycoeff [polyfit $rlist $obslist 2]
   set reslist [polycomp $rlist $obslist $ycoeff]
   set rms [polyrms $reslist]
   echo xx yrms $rms

#Cache xx coeff
   set xxyslope [lindex $ycoeff 1]
   set xxoff [lindex $ycoeff 0]

#Now compute full residuals
   set rms 0.
   set nrms 0
   pgEnv [expr -1.*abs($xrad)] [expr abs($xrad)] -.02 .02 0 0
   update
   foreach x $xlist {
	foreach y $ylist {
	   set xn [expr $x-$xavg]
	   set yn [expr $y-$yavg]
	   set comp [expr $xxoff + $xn*$xxxslope + $yn*$xxyslope]
	   set res [expr $xx($x,$y) - $comp]
	   set rms [expr $rms + pow($res,2)]
	   incr nrms
	   pgPoint $x $res 3
	   }
	}
   set rms [expr sqrt($rms/(1.*$nrms))]
   echo Global rms $rms
   update
   sleep $SLEEP

#Now repeat for xy

#Fit xy in x direction
   set rlist ""
   set obslist ""
   foreach x $xlist {
	lappend rlist [expr $x-$xavg]
	set sum 0.
	foreach y $ylist {
	   set sum [expr $sum + $xy($x,$y)]
	   }
	lappend obslist [expr $sum/(1.*$ny)]
	}
   set xcoeff [polyfit $rlist $obslist 2]
   set reslist [polycomp $rlist $obslist $xcoeff]
   set rms [polyrms $reslist]
   echo xy xrms $rms

#Cache xy coeff
   set xyxslope [lindex $xcoeff 1]

#Fit xy in y direction
   set rlist ""
   set obslist ""
   foreach y $ylist {
	lappend rlist [expr $y-$yavg]
	set sum 0.
	foreach x $xlist {
	   set sum [expr $sum + $xy($x,$y)]
	   }
	lappend obslist [expr $sum/(1.*$nx)]
	}
   set ycoeff [polyfit $rlist $obslist 2]
   set reslist [polycomp $rlist $obslist $ycoeff]
   set rms [polyrms $reslist]
   echo xy yrms $rms

#Cache xy coeff
   set xyyslope [lindex $ycoeff 1]
   set xyoff [lindex $ycoeff 0]

#Now compute full residuals
   set rms 0.
   set nrms 0
   pgEnv [expr -1.*abs($xrad)] [expr abs($xrad)] -.02 .02 0 0
   update
   foreach x $xlist {
	foreach y $ylist {
	   set xn [expr $x-$xavg]
	   set yn [expr $y-$yavg]
	   set comp [expr $xyoff + $xn*$xyxslope + $yn*$xyyslope]
	   set res [expr $xy($x,$y) - $comp]
	   set rms [expr $rms + pow($res,2)]
	   incr nrms
	   pgPoint $x $res 3
	   }
	}
   set rms [expr sqrt($rms/(1.*$nrms))]
   echo Global rms $rms
   update
   sleep $SLEEP

#Now repeat for dilution

#Fit dil in x direction
   set rlist ""
   set obslist ""
   foreach x $xlist {
	lappend rlist [expr $x-$xavg]
	set sum 0.
	foreach y $ylist {
	   set sum [expr $sum + $dil($x,$y)]
	   }
	lappend obslist [expr $sum/(1.*$ny)]
	}
   set xcoeff [polyfit $rlist $obslist 2]
   set reslist [polycomp $rlist $obslist $xcoeff]
   set rms [polyrms $reslist]
   echo dil xrms $rms

#Cache dil coeff
   set dilxslope [lindex $xcoeff 1]

#Fit dil in y direction
   set rlist ""
   set obslist ""
   foreach y $ylist {
	lappend rlist [expr $y-$yavg]
	set sum 0.
	foreach x $xlist {
	   set sum [expr $sum + $dil($x,$y)]
	   }
	lappend obslist [expr $sum/(1.*$nx)]
	}
   set ycoeff [polyfit $rlist $obslist 2]
   set reslist [polycomp $rlist $obslist $ycoeff]
   set rms [polyrms $reslist]
   echo dil yrms $rms

#Cache dil coeff
   set dilyslope [lindex $ycoeff 1]
   set diloff [lindex $ycoeff 0]

#Now compute full residuals
   set rms 0.
   set nrms 0
   update
   pgEnv [expr -1.*abs($xrad)] [expr abs($xrad)] -.02 .02 0 0
   foreach x $xlist {
	foreach y $ylist {
	   set xn [expr $x-$xavg]
	   set yn [expr $y-$yavg]
	   set comp [expr $diloff + $xn*$dilxslope + $yn*$dilyslope]
	   set res [expr $dil($x,$y) - $comp]
	   set rms [expr $rms + pow($res,2)]
	   incr nrms
	   pgPoint $x $res 3
	   }
	}

   set rms [expr sqrt($rms/(1.*$nrms))]

#If noptic = 2, dil is square of dilution, so take sqrt to get arcsec
   if {$noptic == 2} {set rms [expr sqrt($rms)]}
   echo Dilution: rms is [format %.3f $rms] arcsec

#Compute rms whisker error - combined xx and xy.

   set rms2 0.
   set nrms 0
   update
   set whiskmax 0.
   foreach x $xlist {
	foreach y $ylist {
	   set xn [expr $x-$xavg]
	   set yn [expr $y-$yavg]
	   set compxx [expr $xxoff + $xn*$xxxslope + $yn*$xxyslope]
	   set compxy [expr $xyoff + $xn*$xyxslope + $yn*$xyyslope]
	   set resxx [expr $xx($x,$y) - $compxx]
	   set resxy [expr $xy($x,$y) - $compxy]

#Rms whisker error (squared) - note that xy already includes the factor 2.
	   set res2 [expr sqrt(pow($resxx,2) + pow($resxy,2))]
	   set rms2 [expr $rms2 + $res2]

#Whisker itself
	   set whisk [expr sqrt($res2)]
	   set whiskmax [expr max($whiskmax,$whisk)]
	   incr nrms
	   }
	}
   set rms [expr sqrt($rms2/(1.*$nrms))]
   echo Global whisker rms is [format %.3f $rms] arcsec
   echo Max whisker length is [format %.3f $whiskmax] arcsec
   return $rms
   }

######################################################################
#For science report, compute global rms over entire camera.
#args is either 1 or 2 optics designs.

proc wglobal {args} {
   set rms 0.
   set n 0
   foreach ifil "11 12 13 14 18 19 20 24 25 26 30 31 32 35 36 39 40" {
	echo ifil $ifil
	set rms [expr $rms + pow([eval wplot $args $ifil],2)]
	incr n
	}
   set rms [expr sqrt($rms/$n)]
   echo Global rms (whisker length) [format %.3f $rms]
   return $rms
   }

########################################################################
#Check new rtrace diagram parameters.

proc tcheck {hndl xmm ymm ifil} {
   rtrace $hndl $xmm $ymm $ifil 1
   set scale [showFocal $hndl $ifil scale]

#I think this is in mm.
   set fwhmx [exprGet $hndl.diagram->fwhmx]
   set fwhmy [exprGet $hndl.diagram->fwhmy]
   set fwhmxy [exprGet $hndl.diagram->fwhmxy]
   set xx [expr pow($fwhmx*$scale,2)]
   set yy [expr pow($fwhmy*$scale,2)]
   set xy [expr pow($fwhmxy*$scale,2)]

   set e1 [expr $xx - $yy]
   set e2 [expr 2.*$xy]
   set e3 [expr ($xx+$yy)/2.]
   echo rtrace: e1 = $e1 e2 = $e2 e3 = $e3

   set list [polarize $hndl $xmm $ymm $ifil]
   set e1 [lindex $list 0]
   set e2 [lindex $list 1]
   set e3 [lindex $list 2]
   echo polarize: e1 = $e1 e2 = $e2 e3 = $e3
   return
   }

