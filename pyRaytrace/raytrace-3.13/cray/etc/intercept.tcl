#Ray intercept curves

#Make a plot of meridion rays
proc meridianPlot {optic rad filters} {

   set scale 0.
   foreach filter $filters {

#Chief ray
	ray $optic $rad 0 0 0 $filter 0
	set n [exprGet $optic.diagram->np]
	set n1 [expr $n-1]

#Position in focal plane
	set xc [exprGet $optic.diagram->xray<$n1>]
	set yc [exprGet $optic.diagram->yray<$n1>]
	set zc [exprGet $optic.diagram->zray<$n1>]

#Meridional rays
	set xpts($filter) ""
	set ypts($filter) ""
	loop i 0 41 {
	   set xfract [expr ($i-20.)/20.]
	   ray $optic $rad 0 $xfract 0 $filter 0
	   set xf [exprGet $optic.diagram->xray<$n1>]
	   set dx [expr $xf-$xc]
	   lappend xpts($filter) $xfract
	   lappend ypts($filter) $dx
	   set scale [expr max($scale,abs($dx))]
	   }
	}
   set scale [expr 1.1*$scale]
   plotInit a
   pgVport
   pgWindow -1 1 -$scale $scale
   pgBox ant 0 0 ant 0 0
   pgLabel "X Pupil" "DX" "Meridian rays. Radius=$rad filter(s) $filters"
   set sci 1
   set y [expr .9*$scale]
   foreach filter $filters {
	pgLine $xpts($filter) $ypts($filter)
	pgText .9 $y $filter
	set y [expr $y-.1*$scale]
	incr sci
	pgSci $sci
	if {$sci > 7} {set sci 1}
	}
   return
   }

######################################################################
#Make a plot of sagittal rays
proc sagittalPlot {optic rad filters} {

   set scale 0.
   foreach filter $filters {

#Chief ray
	ray $optic $rad 0 0 0 $filter 0
	set n [exprGet $optic.diagram->np]
	set n1 [expr $n-1]

#Position in focal plane
	set xc [exprGet $optic.diagram->xray<$n1>]
	set yc [exprGet $optic.diagram->yray<$n1>]
	set zc [exprGet $optic.diagram->zray<$n1>]

#Sagittal rays.  I will plot dy and dx both (Olso normally does just dy).
	set xpts($filter) ""
	set ypts($filter) ""
	set zpts($filter) ""
	set tpts($filter) ""
	loop i 0 21 {
	   set yfract [expr ($i)/20.]
	   ray $optic $rad 0 0 $yfract $filter 0
	   set xf [exprGet $optic.diagram->xray<$n1>]
	   set dx [expr $xf-$xc]
	   set yf [exprGet $optic.diagram->yray<$n1>]
	   set dy [expr $yf-$yc]
	   lappend xpts($filter) $yfract
	   lappend ypts($filter) $dy
	   lappend zpts($filter) [expr -1.*$yfract]
	   lappend tpts($filter) $dx
	   set scale [expr max($scale,abs($dx))]
	   set scale [expr max($scale,abs($dy))]
	   }
	}
   set scale [expr 1.1*$scale]
   plotInit b
   pgVport
   pgWindow -1 1 -$scale $scale
   pgBox ant 0 0 ant 0 0
   pgLabel "Y Pupil" "DX (Y<0) / DY (Y>0)" \
	"Sagittal rays. Radius=$rad filter(s) $filters"
   set sci 1
   pgSci $sci
   set y [expr .9*$scale]
   foreach filter $filters {
	pgLine $xpts($filter) $ypts($filter)
	pgLine $zpts($filter) $tpts($filter)
	pgText .9 $y $filter
	set y [expr $y-.1*$scale]
	incr sci
	pgSci $sci
	if {$sci > 7} {set sci 1}
	}
   return
   }

#########################################################################
#Make a plot of astigmatism and coma - this is just a recasting of the
#meridian plot above.

proc aberrPlot {optic rad filters} {

   set scale 0.
   foreach filter $filters {

#Chief ray
	ray $optic $rad 0 0 0 $filter 0
	set n [exprGet $optic.diagram->np]
	set n1 [expr $n-1]

#Position in focal plane
	set xc [exprGet $optic.diagram->xray<$n1>]
	set yc [exprGet $optic.diagram->yray<$n1>]
	set zc [exprGet $optic.diagram->zray<$n1>]

#Meridional rays
	set xpts1($filter) ""
	set xpts2($filter) ""
	set astig($filter) ""
	set coma($filter) ""
	loop i 0 21 {
	   set xfract1 [expr ($i)/20.]
	   set xfract2 [expr -1.*$xfract1]
	   ray $optic $rad 0 $xfract1 0 $filter 0
	   set xf1 [exprGet $optic.diagram->xray<$n1>]
	   set dx1 [expr $xf1-$xc]
	   lappend xpts1($filter) $xfract1

	   ray $optic $rad 0 $xfract2 0 $filter 0
	   set xf2 [exprGet $optic.diagram->xray<$n1>]
	   set dx2 [expr $xf2-$xc]
	   lappend xpts2($filter) $xfract2
	   set s [expr ($dx1-$dx2)/2.]
	   set c [expr ($dx1+$dx2)/2.]
	   lappend astig($filter) $s
	   lappend coma($filter) $c
	   set scale [expr max($scale,abs($s))]
	   set scale [expr max($scale,abs($c))]
	   }
	}
   set scale [expr 1.1*$scale]
   plotInit a
   pgVport
   pgWindow -1 1 -$scale $scale
   pgBox ant 0 0 ant 0 0
   pgLabel "X Pupil" "Coma (X<0) / Astimatism (X>0)" \
	"Aberrations. Radius=$rad filter(s) $filters"
   set sci 1
   set y [expr .9*$scale]
   foreach filter $filters {
	pgLine $xpts1($filter) $astig($filter)
	pgLine $xpts2($filter) $coma($filter)
	pgText .9 $y $filter
	set y [expr $y-.1*$scale]
	incr sci
	pgSci $sci
	if {$sci > 7} {set sci 1}
	}
   return
   }

#############################################################################
#Focus.  Compute tangential and sagittal foci vs. field angle
proc focusPlot {optic filters} {

   set scale 0.
   set xfract .02
   foreach filter $filters {
	set xrad [showFocal $optic $filter xrad]
	set field($filter) ""
	set tang($filter) ""
	set sag($filter) ""
	loop i 0 20 {
	   set rad [expr $xrad*$i/20.]

#Chief ray
	   ray $optic $rad 0 0 0 $filter 0
	   set n [exprGet $optic.diagram->np]
	   set n1 [expr $n-1]

#Position in focal plane
	   set xc [exprGet $optic.diagram->xray<$n1>]
	   set yc [exprGet $optic.diagram->yray<$n1>]
	   set zc [exprGet $optic.diagram->zray<$n1>]

#Tangential rays
#Oblique ray
	   ray $optic $rad 0 $xfract 0 $filter 0
	   set n [exprGet $optic.diagram->np]
	   set n1 [expr $n-1]

#Position in focal plane
	   set xf1 [exprGet $optic.diagram->xray<$n1>]
	   set yf1 [exprGet $optic.diagram->yray<$n1>]
	   set zf1 [exprGet $optic.diagram->zray<$n1>]

#Slope at focal plane
	   set mx1 [exprGet $optic.diagram->slope<0>]
	   set my1 [exprGet $optic.diagram->slope<1>]
	   set mz1 [exprGet $optic.diagram->slope<2>]

#Oblique ray other side
	   ray $optic $rad 0 -$xfract 0 $filter 0
	   set n [exprGet $optic.diagram->np]
	   set n1 [expr $n-1]

#Position in focal plane
	   set xf2 [exprGet $optic.diagram->xray<$n1>]
	   set yf2 [exprGet $optic.diagram->yray<$n1>]
	   set zf2 [exprGet $optic.diagram->zray<$n1>]

#Slope at focal plane
	   set mx2 [exprGet $optic.diagram->slope<0>]
	   set my2 [exprGet $optic.diagram->slope<1>]
	   set mz2 [exprGet $optic.diagram->slope<2>]

#Intercept
	   set s1 [expr ($mx2*($zf2 - $zf1) - $mz2*($xf2 - $xf1)) / \
		($mx2*$mz1 - $mx1*$mz2)]
	   set dz [expr $zf1 - $zc + $mz1*$s1]

	   lappend field($filter) $rad
	   lappend tang($filter) $dz
	   set scale [expr max($scale,abs($dz))]

#Now the sagittal rays
#Oblique ray
	   ray $optic $rad 0 0 $xfract $filter 0
	   set n [exprGet $optic.diagram->np]
	   set n1 [expr $n-1]

#Position in focal plane
	   set xf [exprGet $optic.diagram->xray<$n1>]
	   set yf [exprGet $optic.diagram->yray<$n1>]
	   set zf [exprGet $optic.diagram->zray<$n1>]

#Slope at focal plane
	   set mx [exprGet $optic.diagram->slope<0>]
	   set my [exprGet $optic.diagram->slope<1>]
	   set mz [exprGet $optic.diagram->slope<2>]

#y = 0 intercept
	   set s [expr -($yf-$yc)/$my]
	   set dz [expr $zf - $zc + $mz*$s]

	   lappend sag($filter) $dz
	   set scale [expr max($scale,abs($dz))]
	   }
	}

   plotInit b
   pgVport
   pgWindow 0 $rad -$scale $scale
   pgBox ant 0 0 ant 0 0
   pgLabel "Field" "Defocus" "Astigmatic focus filter(s) $filters"
   set sci 1
   set y [expr .9*$scale]
   foreach filter $filters {
	pgSls 1
	pgLine $field($filter) $tang($filter)
	pgSls 2
	pgLine $field($filter) $sag($filter)
	pgSls 2
	pgText [expr .9*$rad] $y $filter
	set y [expr $y-.1*$scale]
	incr sci
	pgSci $sci
	if {$sci > 7} {set sci 1}
	}
   return
   }

