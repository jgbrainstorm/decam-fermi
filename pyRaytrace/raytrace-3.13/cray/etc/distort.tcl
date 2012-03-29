#Compute optical distortions
#Incidence angle stuff has been moved to "tilt.tcl"

#########################################################
#Distortions appropriate for drift scanning

proc distortdrift {hndl filter} {
   global xtracelist xtracelist
   global xofflist xofflist
   global yofflist yofflist
#
   set xtracelist ""
   set xofflist ""
   set yofflist ""
   if {[pgQinf STATE] != "OPEN"} then {pgBegin; pgAsk 0}
#Y limits are 1 pixel
   pgEnv -25.4 25.4 -.024 .024 0 0
   pgLabel "X position (mm)" "Residual (mm)" \
	"Imaging distortion (filter $filter)"
   set xsize [fabs [showparam $hndl -2 [expr $filter+100]]]
   set ysize [fabs [showparam $hndl -1 [expr $filter+100]]]
   for {set i -5} {$i <= 5} {incr i} {
	set xfocal [expr $i*$xsize/5]
	set xoffavg 0
	set yoffavg 0
   for {set j -5} {$j <= 5} {incr j} {
	set yfocal [expr $j*$ysize/5]
	rtrace $hndl $xfocal $yfocal $filter 1
	set xnom [handleShow $hndl.diagram->xnom]
	set ynom [handleShow $hndl.diagram->ynom]
	set xcen [handleShow $hndl.diagram->xcen]
	set ycen [handleShow $hndl.diagram->ycen]
	set xref [handleShow $hndl.diagram->xref]
	set yref [handleShow $hndl.diagram->yref]
#Centroid relative to nominal
	set xoff [expr $xcen-$xnom]
	set yoff [expr $ycen-$ynom]
	set xoffavg [expr $xoffavg+$xoff]
	set yoffavg [expr $yoffavg+$yoff]
	}
	set xoffavg [expr $xoffavg/11]
	set yoffavg [expr $yoffavg/11]
	pgPoint $xfocal $xoffavg 2
	pgPoint $xfocal $yoffavg 4
	lappend xtracelist $xfocal
	lappend xofflist $xoffavg
	lappend yofflist $yoffavg
	}
   }

#########################################################
#Distortions appropriate for spectroscopic fields

proc distort {hndl filter} {
   global xtracelist xtracelist
   global xofflist xofflist
   global yofflist yofflist
#
   set xtracelist ""
   set xofflist ""
   set yofflist ""
   plotInit a
   set xsize [showparam $hndl -2 [expr 100+$filter]]
   set ysize [showparam $hndl -1 [expr 100+$filter]]
   if {$xsize < 0 && $ysize < 0} then {

#Rectangular field
        set rad [sqrt [expr $xsize*$xsize+$ysize*$ysize]]
        } else {

#Circular (or at least elliptical) field
        set rad [max $xsize $ysize]
        }
   pgEnv -$rad $rad -2 2 0 0
   pgLabel "Radius (mm)" "Residual error (mm)" "Spectroscopic distortion \
(filter $filter)"
   for {set i -20} {$i <= 20} {incr i} {
	set xfocal [expr $i*$rad/20]
	set yfocal 0
	rtrace $hndl $xfocal $yfocal $filter 1

#Nominal position is input position.  We convert to angle on sky
#using nominal scale (and 3rd order distortion term, if supplied).
#Xcen and ycen are actual position in focal plane for this angle.
#xoff, yoff are offsets between actual position in focal plane and
#nominal position given by angle on the sky and nominal scale factor.
	set xnom [handleShow $hndl.diagram->xnom]
	set ynom [handleShow $hndl.diagram->ynom]
	set xcen [handleShow $hndl.diagram->xcen]
	set ycen [handleShow $hndl.diagram->ycen]
	set xref [handleShow $hndl.diagram->xref]
	set yref [handleShow $hndl.diagram->yref]

#Centroid relative to nominal
	set xoff [expr $xcen-$xnom]
	set yoff [expr $ycen-$ynom]
	pgPoint $xfocal $xoff 2
	lappend xtracelist $xfocal
	lappend xofflist $xoff
	lappend yofflist $yoff
	}
   }

########################################################################
#
#Fits coefficients to previously computed distortion measurements.
#Run "distort" first.  ndim is actually the number of terms in the fit;
#the max. dimenion is ndim-1.  Typically, ndim is an even number.
#These coefficients are used to convert from angle on the sky to
#position in the focal plane; the equations are
#r = sum[coeff*(ang/scale0)^n] where scale0 is the nominal scale.

proc distortFit {ndim} {
   global xtracelist xtracelist
   global xofflist xofflist
   global yofflist yofflist
   global xcoefflist xcoefflist
   global ycoefflist ycoefflist
   set result ""
   set xcoefflist [polyfit $xtracelist $xofflist $ndim]

#Print out coeffs
   append result "Term     Xcoeff\n"
   loop i 0 [llength $xcoefflist] {
	append result "[format {%2d  %15s} $i [lindex $xcoefflist $i]]\n"
	}
   set xres [polycomp $xtracelist $xofflist $xcoefflist]
   append result "Rms X residuals are [polyrms $xres] mm\n"
   set ycoefflist [polyfit $xtracelist $yofflist $ndim]
   append result "Term     Ycoeff\n"
   loop i 0 [llength $ycoefflist] {
	append result "[format {%2d  %15s} $i [lindex $ycoefflist $i]]\n"
	}
   set yres [polycomp $xtracelist $yofflist $ycoefflist]
   append result " Rms Y residuals are [polyrms $yres] mm\n"
   if {[pgQinf STATE] != "OPEN"} then {pgBegin; pgAsk 0}
#Find radius limits
   set xmin 1.e6
   set xmax -1.e6
   for {set i 0} {$i < [llength $xtracelist]} {incr i} {
	set xmin [min $xmin [lindex $xtracelist $i]]
	set xmax [max $xmax [lindex $xtracelist $i]]
	}
   pgEnv $xmin $xmax -.024 .024 0 0
   pgLabel "Radius (mm)" "Residual (mm)" "Residuals from Fit: nterm = $ndim"
   pgPoint $xtracelist $xres 2
   pgPoint $xtracelist $yres 4
   return $result
   }

#########################################################################
#Compute total distortion and derivative w.r.t. r
proc xcomp {r} {
   global xcoefflist xcoefflist
   set n [llength $xcoefflist]
   set d 0
   set deriv 0
   loop i 0 $n {
	set coeff [lindex $xcoefflist $i]
	set d [expr $d + $coeff * [pow $r $i]]
	if {$i > 0} then {set deriv [expr $deriv + $coeff * $i * [pow $r \
		[expr $i-1]]]}
	}
   return "$d $deriv"
   }

###########################################################################
proc ycomp {r} {
   global ycoefflist ycoefflist
   set n [llength $ycoefflist]
   set d 0
   set deriv 0
   loop i 0 $n {
	set coeff [lindex $ycoefflist $i]
	set d [expr $d + $coeff * [pow $r $i]]
	if {$i > 0} then {set deriv [expr $deriv + $coeff * $i * [pow $i \
		[expr $i-1]]]}
	}
   return "$d $deriv"
   }

proc distortTypeNew {nterm} {
   set type DISTORT$nterm
   set list [typesList]
   if {[lsearch $list $type] >= 0} then {
	echo Type already exists!
	return
	}
   set schema(1) "$type STRUCT [expr $nterm*16+24] 8"
   set schemaElem(1)(1) {ccdrow int 0 NULL 0}
   set schemaElem(1)(2) {ccdcol int 0 NULL 4}
   set schemaElem(1)(3) {xcen float 0 NULL 8}
   set schemaElem(1)(4) {ycen float 0 NULL 12}
   set schemaElem(1)(5) {rayid int 0 NULL 16}
   set schemaElem(1)(6) {wave float 0 NULL 20}
   set schemaElem(1)(7) "xcoeff double 0 $nterm 24"
   set schemaElem(1)(8) "ycoeff double 0 $nterm [expr 24+8*$nterm]"
   schemaDefine schema schemaElem
   unset schema
   unset schemaElem
   return $type
   }

#### Top level procedure to compute distortions in all imager CCDs #####

proc imageDistort {file nterm} {
   set optic [readset $file]
   set type [distortTypeNew $nterm]
#Main CCDs
   set n 5
   set list [chainNew $type]
pgBegin /XWINDOW
pgAsk 0
#Outer loop through columns
   for {set j 1} {$j <= 6} {incr j} {
#Inner loop through rows
	for {set i 5} {$i >= 1} {set i [expr $i-1]} {
	   incr n
	   distortdrift $optic $n
	   distortFit $nterm
	   set hndl [genericNew $type]
	   chainElementAddByPos $list $hndl
	   handleSet $hndl.ccdrow $i
	   handleSet $hndl.ccdcol $j
	   handleSet $hndl.xcen [showparam $optic  14.$n 3]
	   handleSet $hndl.ycen [showparam $optic  14.$n 4]
	   handleSet $hndl.rayid $n
	   handleSet $hndl.wave [showparam $optic -3 [expr $n+100]]
	   for {set k 0} {$k < $nterm} {incr k} {
		handleSet $hndl.xcoeff<$k> [lindex $xcoefflist $k]
		handleSet $hndl.ycoeff<$k> [lindex $ycoefflist $k]
		}
	   handleDel $hndl
	   }
	}
#Next loop through astrometric CCD's
#Row 9
#Loop through columns
   for {set j 1} {$j <= 5} {incr j} {
	   incr n
	   distortdrift $optic $n
	   distortFit $nterm
	   set hndl [genericNew $type]
	   chainElementAddByPos $list $hndl
	   handleSet $hndl.ccdrow 9
	   handleSet $hndl.ccdcol $j
	   handleSet $hndl.xcen [showparam $optic  14.$n 3]
	   handleSet $hndl.ycen [showparam $optic  14.$n 4]
	   handleSet $hndl.rayid $n
	   handleSet $hndl.wave [showparam $optic -3 [expr $n+100]]
	   for {set k 0} {$k < $nterm} {incr k} {
		handleSet $hndl.xcoeff<$k> [lindex $xcoefflist $k]
		handleSet $hndl.ycoeff<$k> [lindex $ycoefflist $k]
		}
	   handleDel $hndl
	}
#Row 6
#Loop through columns
   for {set j 1} {$j <= 5} {incr j} {
	incr n
	   distortdrift $optic $n
	   distortFit $nterm
	   set hndl [genericNew $type]
	   chainElementAddByPos $list $hndl
	   handleSet $hndl.ccdrow 6
	   handleSet $hndl.ccdcol $j
	   handleSet $hndl.xcen [showparam $optic  14.$n 3]
	   handleSet $hndl.ycen [showparam $optic  14.$n 4]
	   handleSet $hndl.rayid $n
	   handleSet $hndl.wave [showparam $optic -3 [expr $n+100]]
	   for {set k 0} {$k < $nterm} {incr k} {
		handleSet $hndl.xcoeff<$k> [lindex $xcoefflist $k]
		handleSet $hndl.ycoeff<$k> [lindex $ycoefflist $k]
		}
	   handleDel $hndl
	}
#Row 8
#Loop through columns
   for {set j 1} {$j <= 6} {incr j} {
	incr n
	   distortdrift $optic $n
	   distortFit $nterm
	   set hndl [genericNew $type]
	   chainElementAddByPos $list $hndl
	   handleSet $hndl.ccdrow 8
	   handleSet $hndl.ccdcol $j
	   handleSet $hndl.xcen [showparam $optic  14.$n 3]
	   handleSet $hndl.ycen [showparam $optic  14.$n 4]
	   handleSet $hndl.rayid $n
	   handleSet $hndl.wave [showparam $optic -3 [expr $n+100]]
	   for {set k 0} {$k < $nterm} {incr k} {
		handleSet $hndl.xcoeff<$k> [lindex $xcoefflist $k]
		handleSet $hndl.ycoeff<$k> [lindex $ycoefflist $k]
		}
	   handleDel $hndl
	}
#Row 7
#Loop through columns
   for {set j 1} {$j <= 6} {incr j} {
	incr n
	   distortdrift $optic $n
	   distortFit $nterm
	   set hndl [genericNew $type]
	   chainElementAddByPos $list $hndl
	   handleSet $hndl.ccdrow 7
	   handleSet $hndl.ccdcol $j
	   handleSet $hndl.xcen [showparam $optic  14.$n 3]
	   handleSet $hndl.ycen [showparam $optic  14.$n 4]
	   handleSet $hndl.rayid $n
	   handleSet $hndl.wave [showparam $optic -3 [expr $n+100]]
	   for {set k 0} {$k < $nterm} {incr k} {
		handleSet $hndl.xcoeff<$k> [lindex $xcoefflist $k]
		handleSet $hndl.ycoeff<$k> [lindex $ycoefflist $k]
		}
	   handleDel $hndl
	}
   set xtbl [schemaTransAuto $type]
   set tbl [schemaToTbl $list $xtbl -schemaName $type]
   hdrInsertLine $tbl.hdr 999 \
"COMMENT Distortion coefficients for design $file"
   hdrInsertLine $tbl.hdr 999 \
"COMMENT CCDROW and CCDCOL give the position of the CCD in the array"
   hdrInsertLine $tbl.hdr 999 \
"COMMENT RAYID gives the ID of the CCD in the ray trace program"
   hdrInsertLine $tbl.hdr 999 \
"COMMENT WAVE gives the wavelength of that CCD in microns"
   hdrInsertLine $tbl.hdr 999 \
"COMMENT X is in the column direction, Y is in the row direction of"
   hdrInsertLine $tbl.hdr 999 \
"COMMENT a CCD frame."
   hdrInsertLine $tbl.hdr 999 \
"COMMENT If x,y is the nominal position of an object integrated across the CCD"
   hdrInsertLine $tbl.hdr 999 \
"COMMENT during a drift scan (units are -25 to 25 mm), then the true position "
   hdrInsertLine $tbl.hdr 999 \
" is given by"
   hdrInsertLine $tbl.hdr 999 \
"COMMENT x(true) = x + xcoeff(0) + xcoeff(1)*x + xcoeff(2)*x^2 + xcoeff(3)*x^3"
   hdrInsertLine $tbl.hdr 999 \
"COMMENT y(true) = y + ycoeff(0) + ycoeff(1)*x + ycoeff(2)*x^2 + ycoeff(3)*x^3"
   echo Writing output to imageCoeff.fit
   fitsWrite $tbl imageCoeff.fit -pdu MINIMAL -binary
   schemaTransDel $xtbl
   distortPrint imageCoeff.fit imageCoeff.txt
   return "\{list $list\} \{table $tbl\}"
   }

########################################################################
#Compute spectroscopic distortions and write to a file.
#Input file name of design so I can document in spec.par

proc specDistort {file filter nterm} {
   global xcoefflist
   set optic [readset $file]
   distort $optic $filter
   distortFit $nterm
#Nominal scale
   set scale [expr abs([showScale $optic $filter])]

echo Writing out parameter file spec.par
   set fid [open spec.par w]
   puts $fid "#Spectroscopic distortions for model $file"
   puts $fid "#If theta is source angle on sky, then nominal radius r (mm)"
   puts $fid "#in the focal plane is"
   puts $fid "#   r = theta/scale"
   puts $fid "#The true radius is given by"
   puts $fid {#   r(true) = r + sum{xcoeff(i)*r^i}}
   puts $fid "#where the sum is for i = 0 to 9"
   puts $fid "#"
   puts $fid "[format {%-10s %15s} scale $scale]         # arcsec/mm"
   puts $fid "[format {%-10s %15s} parity 1]    # Mapping when viewing the sky"
   set wave [showparam $optic -3 [expr 100+$filter]]
   puts $fid "[format {%-10s %15s} wave $wave]   # Design wavelength"
   puts $fid "[format {%-10s %15s} ncoeff [llength $xcoefflist]]\
	# number of coefficients"
   loop i 0 [llength $xcoefflist] {
	puts $fid [format "%-10s %15s" xcoeff($i) [lindex $xcoefflist $i]]
	}
   close $fid

#We no longer compute plate bending distortions.
   return

#Now make a guess at the plate shape for drilling.  We need to bend in
#the reverse direction quite substantially
#Print out the plate bending distortions.
#This is only an approximation.  I assume that the focal plane shape is
# given by
#	z = a2*r^2 + a4*r^4 +a6*r^6
# where r is the cylindrical radius.
#
# There are 2 distortion corrections.  The first is due to the
# finite plate thickness.
# On the bent plate, the concave surface does not compress to 1st order;
# the convex surface stretches.  When a point is on the top convex surface
# that is at radius r, it  maps onto the concave bottom surface
# at radius r - 0.5*w*(dz/dr) where w is the plate
# thickness.  The minus sign has the right sense given the way that our
# z axis is oriented.  The plate thickeness is claimed somewhere to be 3 mm; the
# maximum slope is of order 1 degree; the correction is of order 30 microns
# max, which must be included.  To sufficient accuracy, we have
#	r(mid) = r - .5*w*[2*a2*r + 4*a4*r^3 + 6*a6*r^5]
# Note that since a2 and a4 are negative, we typically have r(mid) > r.
# Next, we unbend the plate.  If s is linear distance along the plate surface,
# then we have
#	s = r + (2/3)*a2^2*r^3 + (8/5)*a2*a4*r^5 + [(12*a2*a6 + 8*a4^2)/7]*r^6
#	      + (24/7)*a2*a6*r^9 + (36/11)*a6^2*r^11
# with r being r(mid).  Thus, s is the "flat"
# position while r is the "distorted" position.  We normally want to use
# this equation in this form: Given the positions of objects on the sky,
# we first compute the cylindrical radius r(true) above and then the radius
# s that it would have in a flattened plate.
echo Writing out parameter file plate1.par
   set fid [open plate1.par w]
   puts $fid "#Plate bending distortions for model $file"
   puts $fid "#Mapping between telescope focal plane and a flat plate"
   puts $fid "#Let r (mm) be the true (cylindrical) radius in the focal plane,"
   puts $fid "#Let s(mm) be the radius in a flat (unbent) plate."
   puts $fid "#Then:"
   puts $fid {#   s = r + sum{xcoeff(i)*r^i}}
   puts $fid "#where the sum is for i = 0 to 11"
   puts $fid "#"
   puts $fid "[format {%-10s %15s} scale 1]"
   puts $fid "[format {%-10s %15s} parity -1]   # Mapping when viewing the plate"
   puts $fid "[format {%-10s %15s} ncoeff 12]	# number of coefficients"
   set a2 [showparam $optic 10 8]
   set a4 [showparam $optic 10 9]
   set a6 [showparam $optic 10 10]
echo a2 = $a2 a4 = $a4 a6 = $a6
   set xcoeff(0) 0
   set xcoeff(1) 0
   set xcoeff(2) 0
   set xcoeff(3) [expr (2./3.)*$a2*$a2]
   set xcoeff(4) 0
   set xcoeff(5) [expr (8./5.)*$a2*$a4]
   set xcoeff(6) 0
   set xcoeff(7) [expr (8.*$a4*$a4 + 12.*$a2*$a6)/7.]
   set xcoeff(8) 0
   set xcoeff(9) [expr 24.*$a4*$a6/9.]
   set xcoeff(10) 0
   set xcoeff(11) [expr 18.*$a6*$a6/11.]
   loop i 0 12 {
   puts $fid [format "%-10s %15s" xcoeff($i) $xcoeff($i)]
	}
   close $fid

#Note: w is a parameter that used to be input, but now I no longer know
#what it is for.
   set w 1
   incidence $optic 3
   distortFit 6
#The coefficients of the fit to incidence angle vs. radius give the
#slope of the surface.  Get the actual coefficients here
   set a2 [expr [lindex $xcoefflist 1]/2.]
   set a4 [expr [lindex $xcoefflist 3]/4.]
   set a6 [expr [lindex $xcoefflist 5]/6.]
#This should have the proper signs and all.
echo Writing out parameter file plate0.par
   set fid [open plate0.par w]
   puts $fid "#Plate bending distortions for model $file"
   puts $fid "#Mapping between drilling plane and a flat plate"
   puts $fid "#Let r (mm) be the true (cylindrical) radius in the focal plane,"
   puts $fid "#Let s(mm) be the radius in a flat (unbent) plate."
   puts $fid "#Then:"
   puts $fid {#   s = r + sum{xcoeff(i)*r^i}}
   puts $fid "#where the sum is for i = 0 to 11"
   puts $fid "#"
   puts $fid "[format {%-10s %15s} scale 1]"
   puts $fid "[format {%-10s %15s} parity -1]   # Mapping when viewing the plate"
   puts $fid "[format {%-10s %15s} ncoeff 12]	# number of coefficients"
   set xcoeff(0) 0
   set xcoeff(1) [expr -2.*$a2*$w]
   set xcoeff(2) 0
   set xcoeff(3) [expr -4.*$a4*$w + (2./3.)*$a2*$a2]
   set xcoeff(4) 0
   set xcoeff(5) [expr -6.*$a6*$w + (8./5.)*$a2*$a4]
   set xcoeff(6) 0
   set xcoeff(7) [expr (8.*$a4*$a4 + 12.*$a2*$a6)/7.]
   set xcoeff(8) 0
   set xcoeff(9) [expr 24.*$a4*$a6/9.]
   set xcoeff(10) 0
   set xcoeff(11) [expr 18.*$a6*$a6/11.]
   loop i 0 12 {
   puts $fid [format "%-10s %15s" xcoeff($i) $xcoeff($i)]
	}
   close $fid
   }

#Print out the linked list

proc distortPrint {infile outfile} {
   if {$outfile == "stdout"} then {set outId stdout} else {
	set outId [open $outfile w]
	}
   set list [fits2Schema $infile DISTORT4]
   chainInit $list iter
   puts $outId "R/C rayid   xc    yc      wave  xcoeff---------------------------------"
   while {1} {
	set hndl [chainNext iter]
	if {$hndl == ""} break
	set line [format " %1d%1d%5d%8.2f%8.2f%6.2f%11.3e%11.3e%11.3e%11.3e" \
	   [handleShow $hndl.ccdrow] \
	   [handleShow $hndl.ccdcol] \
	   [handleShow $hndl.rayid] \
	   [handleShow $hndl.xcen] \
	   [handleShow $hndl.ycen] \
	   [handleShow $hndl.wave] \
	   [handleShow $hndl.xcoeff<0>] \
	   [handleShow $hndl.xcoeff<1>] \
	   [handleShow $hndl.xcoeff<2>] \
	   [handleShow $hndl.xcoeff<3>]]
	puts $outId $line
	}

   chainInit $list iter
   puts $outId "R/C rayid   xc    yc      wave  x=ycoeff---------------------------------"
   while {1} {
	set hndl [chainNext iter]
	if {$hndl == ""} break
	set line [format " %1d%1d%5d%8.2f%8.2f%6.2f%11.3e%11.3e%11.3e%11.3e" \
	   [handleShow $hndl.ccdrow] \
	   [handleShow $hndl.ccdcol] \
	   [handleShow $hndl.rayid] \
	   [handleShow $hndl.xcen] \
	   [handleShow $hndl.ycen] \
	   [handleShow $hndl.wave] \
	   [handleShow $hndl.ycoeff<0>] \
	   [handleShow $hndl.ycoeff<1>] \
	   [handleShow $hndl.ycoeff<2>] \
	   [handleShow $hndl.ycoeff<3>]]
	puts $outId $line
	}
   if {$outId != "stdout"} then {close $outId}
   }
