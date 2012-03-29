#List most recent tracing of a ray
#I can list absolute or lens-o-centric coords
proc rayList {optic {mode ray}} {
   if {$mode != "ray"} {set mode lens}
   msgInit
   set np [exprGet $optic.diagram->np]
   set icolor [exprGet $optic.diagram->icolor]
   set x1 [exprGet $optic.diagram->x${mode}<0>]
   set y1 [exprGet $optic.diagram->y${mode}<0>]
   set z1 [exprGet $optic.diagram->z${mode}<0>]
   set isurf1 [exprGet $optic.diagram->indx<0>]
   set id1 [surfId $optic $isurf1]
   set n1 [showIndex $optic $id1 $icolor]
   set d 0
   set opd -10000
   msgAdd [format "%-6s%12s%12s%12s%12s%12s%12s" \
	   surfid x y z n delta-opl opl]
   loop i 0 $np {
	set x [exprGet $optic.diagram->x${mode}<$i>]
	set y [exprGet $optic.diagram->y${mode}<$i>]
	set z [exprGet $optic.diagram->z${mode}<$i>]
	set delta [expr sqrt(pow($x-$x1,2) + pow($y-$y1,2) + pow($z-$z1,2))]
	set d [expr $delta+$d]
	set opd [expr $opd + $n1*$delta]
	set isurf [exprGet $optic.diagram->indx<$i>]
	set id [surfId $optic $isurf]
	set n [showIndex $optic $id $icolor]
	msgAdd [format "%-6s%12.4f%12.4f%12.4f%12.4f%12.4f%12.4f" \
	   $id $x $y $z $n $delta $d]
	set x1 $x
	set y1 $y
	set z1 $z
	set n1 $n
	}
   msgAdd OPD (mm) = [format %10.5f $opd]
   return [msgGet]
   }

###################################################################
#Plot the ray pattern

proc rayPatternPlot {optic} {
   plotInit b
   pgSci 1
   pgEnv -1.05 1.05 -1.05 1.05 1 0
   pgLabel xfract yfract "Ray Pattern Layout"
   pgSci 6
   set nfract [exprGet $optic.diagram->nfract]
   loop i 0 $nfract {
	set x [exprGet $optic.diagram->xfract<$i>]
	set y [exprGet $optic.diagram->yfract<$i>]
	pgPoint $x $y 3
	}
   pgSci 1
   return
   }

#####################################################################
#Encircle energy for geometric ray trace.  Just list it.
#Allow multiple colors so I can compute effects of lateral chromatic.

# Useful conversion:
# Raytraces: D80 = 1.59 * FWHM where FWHM is 1-D.
# Gaussian: D80 = 1.53 * FWHM
# Gaussian: D80 = 3.59*sigma
# Gaussian: FWHM = 2.35*sigma.
# Diffraction:  FWHM = .88*wave/diam

proc raySum {optic xmm ymm icolors {stopcheck 1}} {
   msgInit
   set xcen 0.
   set ycen 0.
   set n 0.
   foreach icolor $icolors {
	set weight [showFocal $optic $icolor weight]
	rtrace $optic $xmm $ymm $icolor $stopcheck
	set xcen [expr $xcen + [exprGet $optic.diagram->xcen]*$weight]
	set ycen [expr $ycen + [exprGet $optic.diagram->ycen]*$weight]
	set n [expr $n + $weight]
	}
   set xcen [expr $xcen/$n]
   set ycen [expr $ycen/$n]

   set rlist ""
   foreach icolor $icolors {
	set weight [showFocal $optic $icolor weight]
	rtrace $optic $xmm $ymm $icolor $stopcheck
	set nray [exprGet $optic.diagram->nray]
	loop i 0 $nray {
	   set x [exprGet $optic.diagram->xpoint<$i>]
	   set y [exprGet $optic.diagram->ypoint<$i>]
	   set r [expr sqrt(pow($x-$xcen,2) + pow($y-$ycen,2))]
	   lappend rlist $r
	   }
	}
   set np [llength $rlist]
   if {$np == 0} return

#Print out 80% energy
#Scale is from first color.

   set scale [abs [showScale $optic [lindex $icolors 0]]]
   set sum 0.

#The convention is to use diameter
#Uh oh, I really need a weighted radius, but don't have it available.
   foreach r [lsort -real $rlist] {
	set sum [expr $sum + 1./$np]
	if {$sum >= .8} {
	   set rsec [expr $r*$scale]
	   set d [expr $r*2.]
	   set dsec [expr $rsec*2.]
	   if {[verbose]} {
		msgAdd 80% Diameter [format %.4f $d] mm = [format %.2f \
		   $dsec] arcsec
		}
	   if {[verbose]} {echo [msgGet]}
	   return $dsec
	   }
	}
   return
   }

#####################################################################
#PSF shape.
#Let <xx> = expected value of x^2, etc.
#We have
#	e1 = whisklen^2 * cos(2*theta) = <xx> - <yy>
#	e2 = whisklen^2 * sin(2*theta) = 2*<xy>
#	whisklen^2 = sqrt(e1^2 + e2^2)
#	e3 = dilute^2 = (<xx> + <yy>)/2.
#	theta = .5*atan2(e2, e1)
#Jarvis gives e1, e2, dilute all in pixels.

#Compute residual ellipticity for a spot.
#Allow multiple colors so I can compute effects of lateral chromatic.

#####################################################################
#Compute e1, e2, e3.  These are second moments from which one computes
#whisker, dilution amplitudes.
#
#Outputs are correct for FWHM, arcsec.

proc polarize {optic xmm ymm icolors {stopcheck 1}} {

   set xcen 0.
   set ycen 0.
   set n 0.
   set xxsum 0.
   set xysum 0.
   set yysum 0.

   foreach icolor $icolors {
	set weight [showFocal $optic $icolor weight]
	rtrace $optic $xmm $ymm $icolor $stopcheck
	set x [exprGet $optic.diagram->xcen]
	set y [exprGet $optic.diagram->ycen]
	set xx [expr pow([exprGet $optic.diagram->fwhmx]/2.35,2) + $x*$x]

#Fun with signs.
	set fwhmxy [exprGet $optic.diagram->fwhmxy]
	if {$fwhmxy >= 0.} {
	   set sgn 1.
	} else {
	   set sgn -1.
	   }
	set xy [expr $sgn*pow($fwhmxy/2.35,2) + $x*$y]
	set yy [expr pow([exprGet $optic.diagram->fwhmy]/2.35,2) + $y*$y]

#xcen and ycen are weighted properly here.
	set xcen [expr $xcen + $x*$weight]
	set ycen [expr $ycen + $y*$weight]
	set xxsum [expr $xxsum + $xx*$weight]
	set xysum [expr $xysum + $xy*$weight]
	set yysum [expr $yysum + $yy*$weight]
	set n [expr $n + $weight]
	}
   if {$n == 0.} return
   set xcen [expr $xcen/$n]
   set ycen [expr $ycen/$n]

   foreach icolor $icolors {
	set xxsum [expr $xxsum - $xcen*$xcen*$weight]
	set xysum [expr $xysum - $xcen*$ycen*$weight]
	set yysum [expr $yysum - $ycen*$ycen*$weight]
	}

   set xx [expr $xxsum/$n]
   set xy [expr $xysum/$n]
   set yy [expr $yysum/$n]

#This is e1
   set c2t [expr $xx - $yy]

#This is e2
   set s2t [expr 2.*$xy]

#This is the 1-d dilution
   set amp [expr ($xx+$yy)/2.]

#Convert to arcsec^2.  Convert amp to FWHM^2
   set scale [abs [showScale $optic [lindex $icolors 0]]]
   set e1 [expr $c2t*pow($scale*2.35,2)]
   set e2 [expr $s2t*pow($scale*2.35,2)]
   set e3 [expr $amp*pow($scale*2.35,2)]

   return [list $e1 $e2 $e3]
   }

#########################################################################
#Compute whisker length, orientation.  Also project to x, y components.
#Units are arcsec.  Note that amp is sqrt[ecc*ecc/(2.-ecc*ecc)]
#For small ecc, we have ecc = sqrt(2.*eps)

proc whisker {optic xmm ymm icolors {stopcheck 1}} {
   set list [polarize $optic $xmm $ymm $icolors $stopcheck]
   if {[llength $list] == 0} return
   set e1 [lindex $list 0]
   set e2 [lindex $list 1]
   set e3 [lindex $list 2]
   set amp [expr sqrt(sqrt(pow($e1,2) + pow($e2,2)))]
   set theta [expr atan2($e2, $e1)/2.]

#Convert to x and y components, since I will use these later
   set x [expr $amp*cos($theta)]
   set y [expr $amp*sin($theta)]

   set amp [format %.3f $amp]
   set theta [format %.2f [expr $theta*180./3.141593]]
   return [list $x $y $amp $theta]
   }

######################################################################
#Compute dilution.  This is the 1-d average FWHM.

proc dilution {optic xmm ymm icolors {stopcheck 1}} {
   set list [polarize $optic $xmm $ymm $icolors $stopcheck]
   if {[llength $list] == 0} return
   set e1 [lindex $list 0]
   set e2 [lindex $list 1]
   set e3 [lindex $list 2]
   set amp [expr sqrt($e3)]
   set amp [format %.3f $amp]
   return $amp
   }
