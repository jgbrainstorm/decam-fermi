######################################################################
#Run a trace and plot for one or more filters at one spot position
#I believe that there was a bug here that is now fixed.  diagramplot
#plots points relative to xcen, not xref.

proc spotPlot {hndl xmm ymm filters} {
   set first n
   pgSci 1
   set icolor 1
   foreach filter $filters {
	rtrace $hndl $xmm $ymm $filter 1
	if {$first == "n"} then {
	   set xcen [handleShow $hndl.diagram->xcen]
	   set ycen [handleShow $hndl.diagram->ycen]
	   }
	handleSet $hndl.diagram->xcen $xcen
	handleSet $hndl.diagram->ycen $ycen
	handleSet $hndl.diagram->xcen $xcen
	handleSet $hndl.diagram->ycen $ycen

#For multiple filters, space the lines
	if {[llength $filters] > 1} {echo}
	diagramplot $hndl $first
	incr icolor
	pgSci $icolor
	set first y
	}
   if {[llength $filters] > 1} {
	set d80 [raySum $hndl $xmm $ymm $filters]
	set fwhm [expr $d80/1.59]
	echo Combined FWHM = [format %.2f $fwhm]
	}
   }

####################################################################
#Run a grid of spots and make 9 small plots.

proc manyPlot {hndl filters} {
   global XMIN XMAX YMIN YMAX

#Set default to 100 microns, not 1.5 arcsec
   set defaultmm .1
   set scale [expr abs([showFocal $hndl [lindex $filters 0] scale])]
#   set default 1.5
   set default [expr $scale*$defaultmm]
   if {[info exists XMIN]} {set xmin $XMIN} else {set xmin -$default}
   if {[info exists XMAX]} {set xmax $XMAX} else {set xmax $default}
   if {[info exists YMIN]} {set ymin $YMIN} else {set ymin -$default}
   if {[info exists YMAX]} {set ymax $YMAX} else {set ymax $default}
   pgSci 1
   plotInit a 3 3
   set firstpos 0
   foreach ypos "1 0 -1" {
	foreach xpos "-1 0 1" {
	   if {$firstpos == 0} {
		set firstpos 1
	   } else {
		pgNext
		}
	   set icolor 1
	   set first n
	   foreach filter $filters {
		set scale [expr abs([showFocal $hndl $filter scale])]
		set xminmm [expr $xmin/$scale]
		set yminmm [expr $ymin/$scale]
		set xmaxmm [expr $xmax/$scale]
		set ymaxmm [expr $ymax/$scale]
		pgWindow $xminmm $xmaxmm $yminmm $ymaxmm
		set xwid [showFocal $hndl $filter xsize]
		set ywid [showFocal $hndl $filter ysize]
		set xmm [expr $xwid*$xpos]
		set ymm [expr $ywid*$ypos]

#If this is a circle, rescale points so they don't go beyond boundary
		if {$xwid > 0 || $ywid > 0} {
		   set xmm [expr $xmm/1.414]
		   set ymm [expr $ymm/1.414]
		   }
		rtrace $hndl $xmm $ymm $filter 1
		if {$first == "n"} then {
		   set xcen [handleShow $hndl.diagram->xcen]
		   set ycen [handleShow $hndl.diagram->ycen]
		   }
		handleSet $hndl.diagram->xcen $xcen
		handleSet $hndl.diagram->ycen $ycen
		handleSet $hndl.diagram->xcen $xcen
		handleSet $hndl.diagram->ycen $ycen
		pgSci $icolor
		diagramplot $hndl y
		incr icolor
		set first y
		}
	   if {$xpos == 0 && $ypos == 0} {
		pgSci 6
		pgLine "[expr -1.*$xmaxmm/2.] [expr $xmaxmm/2.]" \
		   "$ymaxmm $ymaxmm" 
		set sch [pgQch]
		pgSch 1.5
		pgSci 7
		pgText -[expr $xmaxmm/2.] [expr .8*$ymaxmm] \
		   "[format %.3f $xmaxmm] mm"
		pgText -[expr $xmaxmm/2.] [expr .6*$ymaxmm] \
		   "[format %.1f [expr $xmaxmm*$scale]] arcsec"
		pgSch $sch
		}
	   }
	}
   return
   }

######################################################################
proc spotplot {hndl xmm ymm filters} {
   spotPlot $hndl $xmm $ymm $filters
   return
   }

######################################################################
#Offset a parameter and make a series of spot diagrams.
proc offsetPlot {hndl filter surf param offset} {
   set value [showSurf $hndl $surf $param]
   set newval [expr $value + $offset]
   setSurf $hndl $surf $param $newval
   manyPlot $hndl $filter
   setSurf $hndl $surf $param $value
   return
   }

######################################################################
#Plot the latest spot diagram in an OPTIC structure
#Points are plotted relative to centroid
proc diagramplot {hndl {overlay n}} {
   if {[pgQinf STATE] == "CLOSED" || $overlay == "n"} {
	plotInit a
	}
   set nray [handleShow $hndl.diagram->nray]
   if {$nray == 0} then return
   pgSch 1.0
   global XMIN XMAX YMIN YMAX

#Set default to 100 microns, not 1.5 arcsec
   set defaultmm .1
   set scale [showScale $hndl [exprGet $hndl.diagram->icolor]]
   set scale [expr abs($scale)]
#   set default 1.5
   set default [expr $scale*$defaultmm]

   set xmin -$default; set xmax $default; set ymin -$default;
	set ymax $default

#XMIN, etc are in arcsec.
   if {[info exists XMIN]} {set xmin [expr 1.*$XMIN]}
   if {[info exists XMAX]} {set xmax [expr 1.*$XMAX]}
   if {[info exists YMIN]} {set ymin [expr 1.*$YMIN]}
   if {[info exists YMAX]} {set ymax [expr 1.*$YMAX]}

#Limits in mm
   set xminmm [expr $xmin/$scale]
   set xmaxmm [expr $xmax/$scale]
   set yminmm [expr $ymin/$scale]
   set ymaxmm [expr $ymax/$scale]
   set xcen [handleShow $hndl.diagram->xcen]
   set ycen [handleShow $hndl.diagram->ycen]
   if {$overlay == "n"} then {
	pgPage
	pgSch 1
	pgVport .15 .85 .15 .85
	pgWindow $xmin $xmax $ymin $ymax
	pgSci 5
	pgBox cmst 0 0 cmstv 0 0
	pgSch 1.2
	pgMtext t 3.0 .4 0 "X(arcsec)"
	pgMtext r 3.5 .4 0 "Y(arcsec)"
	pgSch 1
	pgWindow $xminmm $xmaxmm $yminmm $ymaxmm
	pgSci 1
	pgBox bnst 0 0 bnstv 0 0
	pgSch 1.2
	pgMtext b 3.0 .4 0 "X(mm)"
	pgMtext l 3.5 .4 0 "Y(mm)"
	pgSch .8
	pgText [expr $xminmm + .75*($xmaxmm-$xminmm)] \
	   [expr $yminmm+.9*($ymaxmm-$yminmm)] "x = [format %.2f $xcen]"
	pgText [expr $xminmm + .75*($xmaxmm-$xminmm)] \
	   [expr $yminmm+.85*($ymaxmm-$yminmm)] "y = [format %.2f $ycen]"
	}
   set fwhmx [handleShow $hndl.diagram->fwhmx]
   set fwhmy [handleShow $hndl.diagram->fwhmy]
   set fwhmrad [expr sqrt(($fwhmx*$fwhmx + $fwhmy*$fwhmy)/2.)]
   set sigx [expr $fwhmx/2.35]
   set sigy [expr $fwhmy/2.35]
   set sigrad [expr sqrt($sigx*$sigx + $sigy*$sigy)]
   echo [format "x-fwhm(mm) =%12.4f   y-fwhm(mm) =%12.4f   radius = %12.4f" \
	$fwhmx $fwhmy $fwhmrad]
   echo [format "x-sig(mm)  =%12.4f   y-sig(mm)  =%12.4f   radius = %12.4f" \
	$sigx $sigy $sigrad]
   set fwhmx [expr $fwhmx*$scale]
   set fwhmy [expr $fwhmy*$scale]
   set fwhmrad [expr $fwhmrad*$scale]
   set sigx [expr $fwhmx/2.35]
   set sigy [expr $fwhmy/2.35]
   set sigrad [expr sqrt($sigx*$sigx + $sigy*$sigy)]
   echo [format "x-fwhm(arcsec) =%6.2f     y-fwhm(arcsec) =%6.2f     radius =\
%10.2f" $fwhmx $fwhmy $fwhmrad]
   echo [format "x-sig(arcsec)  =%6.2f     y-sig(arcsec)  =%6.2f     radius =\
%10.2f" $sigx $sigy $sigrad]
   set icolor [handleShow $hndl.diagram->icolor]
   pgSch .4
   for {set i 1} {$i < $nray} {incr i} {
	set xpoint [handleShow $hndl.diagram->xpoint<$i>]
	set ypoint [handleShow $hndl.diagram->ypoint<$i>]
	set x [expr ($xpoint-$xcen)]
	set y [expr ($ypoint-$ycen)]
	pgPoint $x $y 3
	update
	}

#Draw airy disk
   set diam [telDiam $hndl]
   set wave [showWave $hndl $icolor]
   set theta [expr 1.22*206265.*$wave*1.e-3/$diam]
   set xlist ""
   set ylist ""
   loop i 0 100 {
	set ang [expr 2.*3.1416*$i/99.]
	set x [expr $theta*cos($ang)/$scale]
	set y [expr $theta*sin($ang)/$scale]
	lappend xlist $x
	lappend ylist $y
	}
   pgSci 6
   pgSls 2
   pgLine $xlist $ylist
   pgSls 1
   pgSci 1
   return
   }

######################################################################
#Plot focal plane layout
proc focalPlot {hndl {colormin 1} {colormax -1}} {
   colorcount $hndl
   set ncolor [handleShow $hndl.ncolor]
   set xmin 1.e6
   set ymin 1.e6
   set xmax 0
   set ymax 0
   set imin [expr $colormin]
   if {$imin < 1} then {set imin 1}
   set imax $colormax
   if {$imax < 0} then {set imax $ncolor}
   if {$imax > $ncolor} then {set imax $ncolor}
   for {set i $imin} {$i <= $imax} {incr i} {
	set flag [showIndex $hndl 0 $i]
	if {$flag == 0} continue
	set scale [expr abs([showFocal $hndl $i scale])]
	set x0 [showFocal $hndl $i xoff]
	set y0 [showFocal $hndl $i yoff]

#Convert arcmin to mm
	if {$scale == 0} {echo Bad scale for i = $i; continue}
	set x0 [expr $x0*60./$scale]
	set y0 [expr $y0*60./$scale]
	set dx [fabs [showFocal $hndl $i xsize]]
	set dy [fabs [showFocal $hndl $i ysize]]

#Plot 4 corners
	set x [expr $x0-$dx]
	set y [expr $y0-$dy]
	set xmax [max $x $xmax]
	set ymax [max $y $ymax]
	set xmin [min $x $xmin]
	set ymin [min $y $ymin]
	set x [expr $x0+$dx]
	set y [expr $y0-$dy]
	set xmax [max $x $xmax]
	set ymax [max $y $ymax]
	set xmin [min $x $xmin]
	set ymin [min $y $ymin]
	set x [expr $x0+$dx]
	set y [expr $y0+$dy]
	set xmax [max $x $xmax]
	set ymax [max $y $ymax]
	set xmin [min $x $xmin]
	set ymin [min $y $ymin]
	set x [expr $x0-$dx]
	set y [expr $y0+$dy]
	set xmax [max $x $xmax]
	set ymax [max $y $ymax]
	set xmin [min $x $xmin]
	set ymin [min $y $ymin]
	}
#   set size [max [fabs $xmin] [fabs $xmax] [fabs $ymin] [fabs $ymax]]
#   set size [expr $size*1.05]
#   set xmin -$size; set xmax $size; set ymin -$size; set ymax $size

#   global XMIN XMAX YMIN YMAX
#   if {[info exists XMIN]} {set xmin $XMIN}
#   if {[info exists XMAX]} {set xmax $XMAX}
#   if {[info exists YMIN]} {set ymin $YMIN}
#   if {[info exists YMAX]} {set ymax $YMAX}
   set xsize [expr $xmax-$xmin]
   set xmin [expr $xmax-1.05*$xsize]
   set xmax [expr $xmin+1.1*$xsize]
   set ysize [expr $ymax-$ymin]
   set ymin [expr $ymax-1.05*$ysize]
   set ymax [expr $ymin+1.1*$ysize]
	
   plotInit b
   pgSci 1
   pgSch 1
   pgEnv $xmin $xmax $ymin $ymax 1 0
   pgLabel "X" "Y" "Focal Plane Layout"
   for {set i $imin} {$i <= $imax} {incr i} {
	set flag [showIndex $hndl 0 $i]
	if {$flag == 0} continue
	set scale [expr abs([showFocal $hndl $i scale])]
	set x0 [showFocal $hndl $i xoff]
	set y0 [showFocal $hndl $i yoff]

#Convert arcmin to mm
	if {$scale == 0} {echo Bad scale for i = $i; continue}
	set x0 [expr $x0*60./$scale]
	set y0 [expr $y0*60./$scale]
	set dx [showFocal $hndl $i xsize]
	set dy [showFocal $hndl $i ysize]

#Plot 4 corners
#	set x [expr $x0-$dx]
#	set y [expr $y0-$dy]
#	pgMove $x $y
#	set x [expr $x0+$dx]
#	set y [expr $y0-$dy]
#	pgDraw $x $y
#	set x [expr $x0+$dx]
#	set y [expr $y0+$dy]
#	pgDraw $x $y
#	set x [expr $x0-$dx]
#	set y [expr $y0+$dy]
#	pgDraw $x $y
#	set x [expr $x0-$dx]
#	set y [expr $y0-$dy]
#	pgDraw $x $y
	pgSfs 2
	if {$dx < 0 || $dy < 0} {
	   pgRect [expr $x0+$dx] [expr $x0-$dx] [expr $y0+$dy] [expr $y0-$dy]
	} else {
	   pgCirc $x0 $y0 $dx
	   }
	pgSch .8
	pgText $x0 $y0 $i .5
	pgSch 1.
	update idletasks
	}
   }

######################################################################
#Compute rms differences in centroids
proc lateral {hndl xmm ymm filters} {
   set xcen 0; set  ycen 0; set nfil 0
   foreach filter $filters {
	rtrace $hndl $xmm $ymm $filter 1
	set xcen [expr $xcen+[handleShow $hndl.diagram->xcen]]
	set ycen [expr $ycen+[handleShow $hndl.diagram->ycen]]
	incr nfil
	}
   set xcen [expr $xcen/$nfil]
   set ycen [expr $ycen/$nfil]
   set xrms 0; set yrms 0
   foreach filter $filters {
	rtrace $hndl $xmm $ymm $filter 1
	set xrms [expr $xrms+[pow [expr \
		[handleShow $hndl.diagram->xcen]-$xcen] 2]]
	set yrms [expr $yrms+[pow [expr \
		[handleShow $hndl.diagram->ycen]-$ycen] 2]]
	}
   set xrms [sqrt [expr $xrms/$nfil]]
   set yrms [sqrt [expr $yrms/$nfil]]
   return "$xrms $yrms"
   }

######################################################################
#### Plot lateral color for a range of wavelengths and positions
proc lateralPlot {hndl filters} {
   set wavemin 9999; set wavemax 0
   foreach filter $filters {
	set wavemin [min $wavemin [showFocal $hndl $filter wave]]
	set wavemax [max $wavemax [showFocal $hndl $filter wave]]
	}
   pgEnv $wavemin $wavemax -.02 .02 0 0
   pgLabel "Wavelength (microns)" "Residual (mm)" \
	"Lateral color for filters $filters"
   set rmax [abs [showFocal $hndl [lindex $filters 0] xsize]]
   loop i 0 5 {
	set xmm [expr $i*$rmax/4.]
	set ymm 0
	set xcen 0; set  ycen 0; set nfil 0;
	foreach filter $filters {
	   rtrace $hndl $xmm $ymm $filter 1
	   set xcen [expr $xcen+[handleShow $hndl.diagram->xcen]]
	   set ycen [expr $ycen+[handleShow $hndl.diagram->ycen]]
	   incr nfil
	   }
	set xcen [expr $xcen/$nfil]
	set ycen [expr $ycen/$nfil]
	set xlist ""; set ylist ""
	foreach filter $filters {
	   rtrace $hndl $xmm $ymm $filter 1
	   set xres [expr [handleShow $hndl.diagram->xcen]-$xcen]
	   set yres [expr [handleShow $hndl.diagram->ycen]-$ycen]
	   lappend xlist [showFocal $hndl $filter wave]
	   lappend ylist $xres
	   }
	pgLine $xlist $ylist
	echo [format "r = %5.0f,   rms = %8.4f" $xmm \
		[lindex [lateral $hndl $xmm $ymm $filters] 0]]
	}
   }

#############################################################
#Compute the rms radius of a spot for one or more filters
proc rrms {hndl xmm ymm filters} {

#First, compute mean rrrms relative to each color's centroid
   set rrms 0; set nfil 0
   foreach filter $filters {
	rtrace $hndl $xmm $ymm $filter 1
	set xrms [expr [handleShow $hndl.diagram->fwhmx]/2.35]
	set yrms [expr [handleShow $hndl.diagram->fwhmy]/2.35]
	set rrms [expr $rrms + $xrms*$xrms + $yrms*$yrms]
	incr nfil
	}
   set rrms [expr $rrms/$nfil]
#Compute rms centroid differences
   set cenrms [lateral $hndl $xmm $ymm $filters]
#Combined
   set xrms [lindex $cenrms 0]
   set yrms [lindex $cenrms 1]
   set rrms [sqrt [expr $rrms + $xrms*$xrms + $yrms*$yrms]]
   return $rrms
   }

#######################################################################
#Run a trace and plot for a drift scan spot

proc driftplot {hndl xmm filter} {
   set i [expr $filter]
   set dy [fabs [handleShow $hndl.fplane<1>->n<$i>]]
#Plot points relative to nominal position
   set first n
   loop j 0 5 {
	set ymm [expr ($j-2.)*$dy/5.]
	rtrace $hndl $xmm $ymm $filter 1
	handleSet $hndl.diagram->xcen [handleShow $hndl.diagram->xnom]
	handleSet $hndl.diagram->ycen [handleShow $hndl.diagram->ynom]
	diagramplot $hndl $first
	set first y
	}
   }

######################################################################
#Compute and plot the incidence angles for a set of spots in one filter
proc incplot {hndl filter} {
   set xmax [showparam $hndl -2 [expr 100+$filter]]
   set ymax [showparam $hndl -1 [expr 100+$filter]]
   if {$xmax < 0 && $ymax < 0} then {
#Rectangular field
	set rmax [sqrt [expr $xmax*$xmax+$ymax*$ymax]]
	} else {
#Circular (or at least elliptical) field
	set rmax [max $xmax $ymax]
	}
   if {[pgQinf STATE] == "CLOSED"} then {pgBegin; pgAsk 0}
   pgEnv 0 $rmax 0 2. 0 0
   loop i 0 25 {
	set rad [expr $i*$rmax/24.]
	rtrace $hndl $rad 0 $filter 1
	set ang [expr [handleShow $hndl.diagram->incidence]*57.3]
	pgPoint $rad $ang 3
	}
#Alternative computation should get same answer.
   loop i 0 25 {
	set rad [expr $i*$rmax/24.]
	rtrace $hndl $rad 0 $filter 1
	set slope [handleShow $hndl.diagram->slope]
	set normal [handleShow $hndl.diagram->normal]
	set ang [expr ($slope-$normal)*57.3]
	pgPoint $rad $ang 4
	}
   }

#Function diagramSim has been deleted - use spotMap in geom.tcl instead.
#(diagramSim called an astrotools function which wouldn't work anyway!)
