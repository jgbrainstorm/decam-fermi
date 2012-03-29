#Support routines for cheapboss design.

#I have duplicate filters, 1-4 and 5-8.  I will use 5-8 to put in
#refracted positions.

proc cheapRefract {hndl zenith} {

#I will shift the xoff focal plane parameter

#Reference wavelength for refraction
   set wave0 [showFocal $hndl 1 wave]
   foreach ifil "5 6 7 8" {
	set wave [showFocal $hndl $ifil wave]

#Refraction is in arcsec.  xoff is in arcmin
	set refract [diffRefract $wave $wave0 $zenith]
	set scale [showScale $hndl $ifil]
	set offset [expr $refract/$scale]
	set surf [format %02f [expr 15. + $ifil/100.]]
	setSurf $hndl $surf x [format %.3f $offset]
echo surf $surf offset $offset
	}
   return
   }

##################################################################
#Define spot positions and weights
proc despecSpots {} {
   clearSpot
   addSpot 0 0 .2
   addSpot .2 0 .4
   addSpot .40 .6
   addSpot .6 0 .8
   addSpot .8 0 1.
   addSpot 1. 0 1.2
   return
   }

########################################################################
#Define default flags.  Need to add curvature flags on top of these.

proc despecFlags {hndl} {
   clearFlags
   linkColorFlag $hndl "1 2 3 4" 2
   setSurfFlag $hndl 2 z 1
   linkSurfFlag $hndl "17 18" z 4
   setSurfFlag $hndl 19 z 1
   return
   }

########################################################################
#Measure thickness of ADC.  Note that zsurf gives lens-o-centric z position
#and is not useful for surfaces that are tilted.

#This procedure is specific to despec and thus has ADC surface IDs
#hardwired.

proc adcThick {hndl} {
   set lenses "9,10 10,11 12,13 13,14"
   foreach x "-225 0 225" {
	if {$x == 0} {
	   set xfract 0
	} elseif {$x > 0} {
	   set xfract 1
	} else {
	   set xfract -1
	   }
	foreach lens $lenses {
	   set dzmin($x,$lens) 999.
	   set dzmax($x,$lens) -999.
	   }
	ray $hndl $x 0 $xfract 0 1 0
	foreach i "9 10 11 12 13 14 15" {
	   set z($i) [exprGet $hndl.diagram->zray<$i>]
	   }
	foreach lens $lenses {
	   set list [split $lens ,]
	   set l1 [lindex $list 0]
	   set l2 [lindex $list 1]
	   set dz($x,$lens) [expr abs($z($l1) - $z($l2))]
	   }
	}

   echo [format "%-5s  %5s  %8s" Lens X Thick]
   foreach lens $lenses {
	echo ""
	foreach x "-225 0 225" {
	   echo [format "%-5s  %5.0f  %8.3f" $lens $x $dz($x,$lens)]
	   }
	}

#Get total thickness of glued ADC
   set thick1left [expr $dz(-225,9,10) + $dz(-225,10,11)]
   set thick1center [expr $dz(0,9,10) + $dz(0,10,11)]
   set thick1right [expr $dz(225,9,10) + $dz(225,10,11)]
   set thick1 [max $thick1left $thick1center $thick1right]

   set thick2left [expr $dz(-225,12,13) + $dz(-225,13,14)]
   set thick2center [expr $dz(0,12,13) + $dz(0,13,14)]
   set thick2right [expr $dz(225,12,13) + $dz(225,13,14)]
   set thick2 [max $thick2left $thick2center $thick2right]

   echo ADC1 combined thickness [format %.3f $thick1]
   echo ADC2 combined thickness [format %.3f $thick2]

   return
   }

##########################################################################
#Get thinnest thickness on a lens.  Inputs are surface indexes, not names.

proc thinnest {hndl isurf1 isurf2} {
   set thin 999.
   for {set j -20} {$j <= 20} {incr j} {
	set x [expr $j*225./20.]
	for {set i -20} {$i <= 20} {incr i} {
	   set xfract [expr $i*1./20.]
	   ray $hndl $x 0. $xfract 0. 1 0
	   set z1 [exprGet $hndl.diagram->zray<$isurf1>]
	   set z2 [exprGet $hndl.diagram->zray<$isurf2>]
	   set dz [expr abs($z2-$z1)]
	   set thin [expr min($thin,$dz)]
	   }
	}
   return [format %.3f $thin]
   }

########################################################################

#The following stuff runs in kentools
#Create a region with an exponential profile.
#re is in pixels.  I will not attach units (arcsec) yet.
#tot is arbitrary normalization.

proc expSim {re {tot 1000.}} {

#Convert re to h
   set h [expr $re/1.68]
   set nrow [expr round(2.*10.*$h)]
   set nrow [expr (($nrow+1)/2)*2 + 1]
   set reg [regNew -type FL32 $nrow $nrow]
   regSetWithDbl $reg 0.

   if {$h == 0.} {
	regPixSetWithDbl $reg 0 0 $tot
	return
	}

   set sum 0.
   loop i 0 $nrow {
	loop j 0 $nrow {
	   set x [expr $i - ($nrow-1)/2]
	   set y [expr $j - ($nrow-1)/2]
	   set r [expr sqrt($x*$x + $y*$y)]
	   set pix [expr exp(-$r/$h)]
	   regPixSetWithDbl $reg $i $j $pix
	   set sum [expr $sum + $pix]
	   }
	}
   set norm [expr $tot/$sum]
   loop i 0 $nrow {
	loop j 0 $nrow {
	   set pix [expr $norm * [regPixGet $reg $i $j]]
	   regPixSetWithDbl $reg $i $j $pix
	   }
	}
   return $reg
   }

########################################################################
#Create a region with a Gaussian PSF.  I will input FWHM.

proc gaussSim {fwhm} {

#Total intensity is 1000.  Make region go out to 5*h and be odd.
   set tot 1.
   set sig [expr $fwhm/2.35]
   set nrow [expr round(2.*3.*$fwhm)]
   set nrow [expr (($nrow+1)/2)*2 + 1]
   set reg [regNew -type FL32 $nrow $nrow]
   regSetWithDbl $reg 0.

   if {$sig == 0.} {
	regPixSetWithDbl $reg 0 0 $tot
	return
	}

   set sum 0.
   loop i 0 $nrow {
	loop j 0 $nrow {
	   set x [expr $i - ($nrow-1)/2]
	   set y [expr $j - ($nrow-1)/2]
	   set r [expr sqrt($x*$x + $y*$y)]
	   set pix [expr exp(-$r*$r/(2.*$sig*$sig))]
	   regPixSetWithDbl $reg $i $j $pix
	   set sum [expr $sum + $pix]
	   }
	}
   set norm [expr $tot/$sum]
   loop i 0 $nrow {
	loop j 0 $nrow {
	   set pix [expr $norm * [regPixGet $reg $i $j]]
	   regPixSetWithDbl $reg $i $j $pix
	   }
	}
   return $reg
   }

######################################################################
#Convolve two regions - return in a new region.

proc convolve {reg1 reg2} {
   set nrow [exprGet $reg1.nrow]
   set ncol [exprGet $reg1.ncol]
   set regout [regNew -type FL32 $nrow $ncol]
   regSetWithDbl $regout 0.
   regConvolve $reg1 $reg2 -regtarget $regout
   return $regout
   }

########################################################################
#Conclusions regarding 1/2-light radius of COSMOS galaxies:
#
# Radius is log-normal
#       sigma = .2
#       median = 3.66 - .114*mag(ACS)
# Units are log10(ACS pixels) (.03 arcsec)

#Loop through a set of magnitudes, radii, and apertures

proc despecSim {} {
   global img specmag

   if {[info exists specmag]} {unset specmag}

   set img(sky) 0.

#Scale (arcsec/pixel)
   set scale 0.1
   set norm 1000.

   set specmag() "mag rfact fwhm apdiam"
   set specmag(mag) "18 19 20 21 22"
   set specmag(rfact) ".63 1. 1.58"
   set specmag(fwhm) "0.8 1.1 1.4"
   set specmag(apdiam) "1.0 1.5 2.0 2.5 3.0"
   foreach mag $specmag(mag) {

#Median 1/2-light radius in arcsec.  0.03 is size of ACS pixel.
	set median [expr 0.03*pow(10.,3.66 - 0.114*$mag)]
	set medpix [expr $median/$scale]

#Plus and minus 1-sigma of median are factors 1.58, .63
	foreach rfact $specmag(rfact) {
	   set rad [expr $median*$rfact]
	   set radpix [expr $rad/$scale]
	   set reg1 [expSim $radpix $norm]

#Range of FWHM, from best of the best to random.
	   foreach fwhm $specmag(fwhm) {
		set fwhmpix [expr $fwhm/$scale]
		set reg2 [gaussSim $fwhmpix]
		set regout [convolve $reg1 $reg2]
		set center [expr [exprGet $regout.nrow]/2.]
		axSet row $center
		axSet col $center

#Fiber diameter
		foreach apdiam $specmag(apdiam) {
		   set apdiampix [expr $apdiam/$scale]
		   set apradpix [expr $apdiampix/2.]
		   circle $regout $apradpix

#Fraction of light in aperture
		   set cnts $img(phototot)
		   set fract [expr $cnts/$norm]
		   set dmag [expr -2.5*log10($fract)]
		   set apmag [expr $mag + $dmag]
		   set specmag($mag,$rfact,$fwhm,$apdiam) $apmag
		   }
		regDel $regout
		regDel $reg2
		}
	   regDel $reg1
	   }
	}
   return
   }

################################################################

#Plot out results of despecSim
#The main question is choice of aperture.

proc specmagPlot {} {
   global specmag

#i-band sky brightness at secz = 1.3
#Zenith is 20.2, slope is 0.6
   set skymag [expr 20.2 - (1.3-1)*0.6]

   foreach mag $specmag(mag) {
	foreach apdiam $specmag(apdiam) {
	   set medlist($mag,$apdiam) ""
	   }
	}
   foreach mag $specmag(mag) {
	foreach rfact $specmag(rfact) {
	   foreach fwhm $specmag(fwhm) {
		foreach apdiam $specmag(apdiam) {
set ap [expr $apdiam + ([random 100]-99./2.)/1000.]
set apmag $specmag($mag,$rfact,$fwhm,$apdiam)

#Signal/Noise
#Let's do a better job of S/N
		   set sig [expr pow(10.,-$apmag/2.5)]
		   set area [expr 3.14*pow($apdiam,2)/4.]
		   set sky [expr $area*pow(10.,-$skymag/2.5)]

#Normalize sky to signal from total mag of galaxy.
#By doing so, my S/N is the same as the mag. of a galaxy for which we
#collect all the flux and have zero sky.
		   set norm [expr pow(10.,-$mag/2.5)]
		   set sn [expr -2.5*log10($sig/sqrt(($sig+$sky)/$norm))]
		   lappend medlist($mag,$apdiam) $sn
		   }
		}
	   }
	}

#Compute median S/N
   pgEnv 0 3.5 18. 25. 0 0
   pgLabel "Aperture (arcsec)" "Effective Magnitude(small is good)" \
	"Simulated spectroscopy Signal/Noise v. Aperture"
   foreach mag $specmag(mag) {
	set xlist ""
	set ylist ""
	foreach apdiam $specmag(apdiam) {
	   lappend xlist $apdiam
	   set list [lsort -real $medlist($mag,$apdiam)]
	   set median [lindex $list [expr [llength $list]/2]]
	   lappend ylist $median
	   echo Mag $mag Aperture Diameter $apdiam \
		Median S/N [format %.2f $median]
	   }
	pgLine $xlist $ylist
	pgText 3.2 [lindex $ylist end] $mag
	}
   return
   }

################################################################
#Make a set of plots in .gif files of the ADC elements - to be printed.

proc adcPlot {hndl} {
   global XMIN YMIN XMAX YMAX

   if {[info exists XMIN]} {unset XMIN}
   if {[info exists XMAX]} {unset XMAX}
   set YMIN -10450
   set YMAX -10000
   set YL -10400

   pgEnd
   global INVERT
   set INVERT 1
   set SIZE 1200

#ADC1-1
   pgBegin adc1-1.gif/GIF

#Make big - to cut down on jaggies
   pgGeomSet $SIZE $SIZE 0 0
   pgSci 1
   opticPlot h0 9 10 1 no
   pgSci 1
   pgText 0 $YL "ADC1-1"
   pgEnd

#ADC1-2
   pgBegin adc1-2.gif/GIF

#Make big - to cut down on jaggies
   pgGeomSet $SIZE $SIZE 0 0
   pgSci 1
   opticPlot h0 10 11 1 no
   pgSci 1
   pgText 0 $YL "ADC1-2"
   pgEnd

#ADC2-1
   pgBegin adc2-1.gif/GIF

#Make big - to cut down on jaggies
   pgGeomSet $SIZE $SIZE 0 0
   pgSci 1
   opticPlot h0 12 13 1 no
   pgSci 1
   pgText 0 $YL "ADC2-1"
   pgEnd

#ADC2-2
   pgBegin adc2-2.gif/GIF

#Make big - to cut down on jaggies
   pgGeomSet $SIZE $SIZE 0 0
   pgSci 1
   opticPlot h0 13 14 1 no
   pgSci 1
   pgText 0 $YL "ADC2-2"
   pgEnd

   unset YMIN YMAX
   unset INVERT
   return
   }

######################################################################
#Routine to take cuts in surface height across C4 to verify that zernikes
#are doing the right thing.

proc zcut {old new y} {
   set rmax 287
   set isurf [surfIndex $old 11]
   set xlist ""; set zlist ""
   for {set x -$rmax} {$x <= $rmax} {incr x} {
	set r [expr sqrt($x*$x + $y*$y)]
	if {$r > $rmax} continue
	set zold [zsurf $old $isurf $x $y]
	set znew [zsurf $new $isurf $x $y]

#microns
	set zdiff [expr 1.e3*($znew - $zold)]
	lappend xlist $x
	lappend zlist $zdiff
	}
   plotInit a
   pgEnv -$rmax $rmax -5. 5. 0 0
   pgLine $xlist $zlist
   pgLabel "x (mm)" "zdiff (microns)" "Cut at Y = $y"
   return
   }

#########################################################
#Convert r, theta to x,y

proc rt-xy {r theta} {
   set x [expr $r*cos($theta/57.3)]
   set y [expr $r*sin($theta/57.3)]
   return [list $x $y]
   }

#########################################################

#Make a fits image of psf increments due to C4

proc psfInc {old new} {
   set NPIX 199
   set reg [regNew $NPIX $NPIX]
   loop r 0 $NPIX {
	echo row $r out of $NPIX
	loop c 0 $NPIX {
	   set yfract [expr ($r+0.5 - $NPIX/2.)/($NPIX/2.)]
	   set xfract [expr ($c+0.5 - $NPIX/2.)/($NPIX/2.)]
	   set rfract [expr sqrt($xfract*$xfract + $yfract*$yfract)]
	   if {$rfract >= 1.} continue
	   set xmm [expr $xfract*225.]
	   set ymm [expr $yfract*225.]
	   rtrace $old $xmm $ymm 3 1
	   rtrace $new $xmm $ymm 3 1
	   set oldfwhmx [exprGet $old.diagram->fwhmx]
	   set oldfwhmy [exprGet $old.diagram->fwhmy]
	   set newfwhmx [exprGet $new.diagram->fwhmx]
	   set newfwhmy [exprGet $new.diagram->fwhmy]

	   set oldblurx [exprGet $old.diagram->blurx]
	   set oldblury [exprGet $old.diagram->blury]
	   set newblurx [exprGet $new.diagram->blurx]
	   set newblury [exprGet $new.diagram->blury]

	   set fwhmdiff [expr pow($newfwhmx,2) + pow($newfwhmy,2) \
		- pow($oldfwhmx,2) - pow($oldfwhmy,2)]

	   set sign 1.
	   if {$fwhmdiff < 0.} {
		set sign -1.
		set fwhmdiff [expr -1.*$fwhmdiff]
		}

#Convert to microns
	   set fwhmdiff [expr sqrt($fwhmdiff)*1.e3*$sign]

#This is fwhmx + fwhmy.  Divide by 2.35 to get rms radius.
	   set fwhmdiff [expr $fwhmdiff/2.35]

	   regPixSet $reg $r $c $fwhmdiff
	   }
	}
   regWriteAsFits $reg diff.fit
   regDel $reg
   return diff.fit
   }

###################################################################
#Make a plot of the WFE around the rim of a zernike structure.
#I do not have ds9 installed on n150 at the moment ...

proc zrim {zern {rad .99}} {
   set nstep 300
   set tlist ""
   set errlist ""
   loop i 0 $nstep {
   	set theta [expr 360.*($i*1./$nstep)]
	set list [rt-xy $rad $theta]
	set xfract [lindex $list 0]
	set yfract [lindex $list 1]
	set err [zernikeWaveErr $zern $xfract $yfract]
	lappend tlist $theta
	lappend errlist $err
	}
   pgEnv 0 360 -.01 .01 0 0
   pgLine $tlist $errlist
   return
   }

###############################################################
#Take steps in radius along x axis.  Only the cn terms (cn = 0) are relevant
#
proc kent {hndl} {

#We go around the focal plane at constant radius and compute Zernikes
#This routine caches the Zernike terms for later plotting.

#Possible list of m values
   set mlist(2) "0 2"
   set mlist(3) "1 3"
   set mlist(4) "0 2 4"
   global zcoeff
   if {[info exists zcoeff]} {unset zcoeff}
   set zcoeff() ""
   set ang 0.
   loop i 0 365 10 {

#Radius in degrees.
	set rad [expr ($i/360.)*1.1]
	set radmm [expr $rad*3600./17.57]
	lappend zcoeff() $radmm
	set xrad [expr $radmm * cos($ang/57.2958)]
	set yrad [expr $radmm * sin($ang/57.2958)]

#Filter 5 of blanco-2602 is the nearest to Roodman's .7 microns
	set zern [zernikeWaveFit $hndl $xrad $yrad 5]
	set norder [exprGet $zern.norder]
	foreach n "2 3 4" {
	   foreach m $mlist($n) {
		if {$m == 0} {
		   set clist 0
		} else {
		   set clist "0 1"
		   }
		foreach cn $clist {
		   set zcoeff($radmm,$n,$m,$cn) 0.
		   if {$n+$m > $norder} continue
		   set zcoeff($radmm,$n,$m,$cn) [zzshow $zern $n $m $cn]
		   }
		}
	   }
	zernikeDel $zern
	}
   return
   }

###############################################################
#Make plots that match Aaron's Zemax plots.  rad is in degrees
#
proc rood {hndl rad} {

#We go around the focal plane at constant radius and compute Zernikes
#This routine caches the Zernike terms for later plotting.

#Possible list of m values
   set mlist(2) "0 2"
   set mlist(3) "1 3"
   set mlist(4) "0 2 4"
   global zcoeff
   if {[info exists zcoeff]} {unset zcoeff}
   set zcoeff() ""
   loop ang 0 360 15 {
	lappend zcoeff() $ang
	set xrad [expr $rad*3600./17.57 * cos($ang/57.2958)]
	set yrad [expr $rad*3600./17.57 * sin($ang/57.2958)]

#Filter 5 of blanco-2602 is the nearest to Roodman's .7 microns
	set zern [zernikeWaveFit $hndl $xrad $yrad 5]
	foreach n "2 3 4" {
	   foreach m $mlist($n) {
		if {$m == 0} {
		   set clist 0
		} else {
		   set clist "0 1"
		   }
		foreach cn $clist {
		   set zcoeff($ang,$n,$m,$cn) [zzshow $zern $n $m $cn]
		   }
		}
	   }
	zernikeDel $zern
	}
   return
   }

###################################################################
#We've run "rood".

proc roodPlot {n m cn} {
   global zcoeff
   set xlist ""
   set ylist ""
   foreach ang $zcoeff() {
	lappend xlist $ang
	lappend ylist $zcoeff($ang,$n,$m,$cn)
	}
#   lappend xlist [expr [lindex $xlist 0] + 360]
#   lappend ylist [lindex $ylist 0]

   set xmin [eval min $xlist]
   set xmax [eval max $xlist]
   set diff [expr $xmax-$xmin]
   set xmin [expr $xmin - 0.1*$diff]
   set xmax [expr $xmax +0.1*$diff]

   set ymin [eval min $ylist]
   set ymax [eval max $ylist]
   set diff [expr $ymax-$ymin]
   set ymin [expr $ymin - 0.1*$diff]
   set ymax [expr $ymax +0.1*$diff]
   pgEnv $xmin $xmax $ymin $ymax 0 0
   pgLine $xlist $ylist
   pgLabel Angle/Radius "Zerike Coeff Amplitude (waves)" \
	"Zernike Term   n = $n   m = $n"
   return
   }

###########################################################################
#Trim a wavefront error map.

proc regTrim {reg row0 row1 col0 col1} {

#require row1-row0 = col1-col0 (same radius).
   if {$row1-$row0 != $col1-$col0} {
	echo Delta-row and Delta-col must be the same!
	return
	}
   set nrow [expr $row1-$row0+1]
   set ncol [expr $col1-$col0+1]
   set subreg [subRegNew $reg $nrow $ncol $row0 $col0]
   set outreg [regNewFromReg $subreg]
   regDel $subreg

#Now zero out pixels outside the radius.
   set rcen [expr $nrow/2.]
   set ccen [expr $ncol/2.]
   set rmax $rcen
   loop i 0 $nrow {
	loop j 0 $ncol {
	   set rad2 [expr pow(($i+0.5)-$rcen,2) + pow(($j+0.5)-$ccen,2)]
	   if {$rad2 > $rmax*$rmax+1.} {regPixSetWithDbl $outreg $i $j 0.}
	   }
	}
   return $outreg
   }

#######################################################################
#Get a global PSF
#Return global rms FWHM in arcsec

proc psfReport {hndl} {
   set fracts "0 .33 .66 1.0"
   set weights ".1 .4 .7 1."

   set wsum 0.
   set sum 0.
   set filters "1 2 3 4"
   loop i 0 [llength $fracts] {
        set fract [lindex $fracts $i]
        set weight [lindex $weights $i]
        set rmm [expr 225.*$fract]

#raySum returns D80.
        set fwhm [expr [raySum $hndl $rmm 0. $filters]/1.59]
        set sum [expr $sum + $fwhm*$fwhm*$weight]
        set wsum [expr $wsum + $weight]
        }
   set fwhm [expr sqrt($sum/$wsum)]
   return [format %.2f $fwhm]
   }
