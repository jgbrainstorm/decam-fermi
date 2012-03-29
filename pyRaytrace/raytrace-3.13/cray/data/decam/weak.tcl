#Compute spot sizes, including effects of refraction.
#If filter and lens names are set:
#
#filterSample <optic> xmm ymm ref-color series filter
#prismSet <optic1> zenith
#prismSolve <optic1>
#prismReset <optic1>
#
#Top-level:  adcEval <optic> rad ref-color series filter
#
#Whisker analysis:
#	whiskerList "hndls"  Give handles to may11, may11a, adc5, vista
#				Compute diff. refraction effects.
#	whiskLin "hndls"     Give handles to may11, may11a, adc5, vista
#				Compute psf change across focal plane
#	whiskFocus "hndls"   Give handles to designs with steps in focus or
#					other parameters.
#				Compute psf changes between multiple designs.
#
#	whiskPlot h1 h2 ifil	Plot differences in PSF between 2 designs.
#	m2			Plot mosaic II outline
#
#First, set glass types in may11 and adc designs.  Also, add 4 extra colors
#for scratch space.

###################################################################
#Perform a dense sampling of a given filter.

proc denseSample {optic icolor1 icolor2} {

   set wave1 [showWave $optic $icolor1]
   set wave2 [showWave $optic $icolor2]

   set ncolor [exprGet $optic.ncolor]

#Clear out old extra colors.
   if {$ncolor > 8} {
	for {set i 9} {$i <=$ncolor} {incr i} {setIndex $optic 0 $i 0}
	}

   colorcount $optic 

#I don't want full weight at wavelength limits, so I will deposit 4
#wavelengths uniformly covering the interval but staying inside the wavelength
#limits.
   loop i 0 4 {
	set wave [expr $wave1 + ($wave2-$wave1)*($i+.5)/4.]
	waveAdd $optic $wave $icolor1
	}
   return
   }

###################################################################
#Perform a dense sampling of a given filter.
#I read filter response curve from a file.
#To make this general, I will make a copy of the design and implement just
#the wavelengths from the filter file.
#I will use "icolor" as a template for configuration information.

proc filterSample {inoptic icolor series filter} {
   global env

#Find the filter first.  Because I sort filters by series, I can't use
#dataFileFind, which would otherwise work fine if I just supplied the
#name of the file without path info.
   set file $env(CRAY_DIR)/data/filter/$series/$filter.fil
   if {![file exists $file]} {
	error "Unable to find file $series/$filter.fil"
	}
   set waves ""
   set weights ""
   set fid [open $file]
   while {1} {
	set line [string trim [gets $fid]]
	if {[eof $fid]} break
	if {[string length $line] == 0} continue
	if {[string index $line 0] == "#"} continue
	set wave [lindex $line 0]
	set weight [lindex $line 1]
	lappend waves $wave
	lappend weights $weight
	}
   close $fid

   set optic [opticNew]
   opticCopy $inoptic $optic

#Copy template to filter 1
   if {$icolor != 1} {
	waveCopy $optic $icolor 1
	}
   set ncolor [exprGet $optic.ncolor]

#Clear out old extra colors.
   if {$ncolor > 1} {
	for {set i 2} {$i <=$ncolor} {incr i} {setIndex $optic 0 $i 0}
	}

   colorcount $optic 

   set wave [lindex $waves 0]
   set weight [lindex $weights 0]
   waveSwitch $optic 1 $wave
   setFocal $optic 1 weight $weight
   loop i 1 [llength $waves] {
	set wave [lindex $waves $i]
	set weight [lindex $weights $i]
	set icolor [expr $i+1]
	waveAdd $optic $wave 1
	setFocal $optic $icolor weight $weight
	}
   return $optic
   }

#######################################################################
#Set for optimal rotation of prisms.  Assume that I have run filterSample
#This works even for non-adc design although phi angles are meaningless.
#CHANGE - just set yoff.  I used to set prism rotation to nominal values,
#which is needed when I solve for optimal rotation (otherwise the
#centered derivatives are zero!) but I lost track of that fact and was
#assuming my spot diagrams were reflecting uncorrected refraction.  Oops!
#Instead, I set prisms to a nominal offset in prismSolve

proc prismSet {optic zenith} {

   set ncolor [exprGet $optic.ncolor]
   set nsurf [exprGet $optic.ncolor]

#Now set yoff to exactly match refraction
   set wave1 [showWave $optic 1]
   for {set i 2} {$i <= $ncolor} {incr i} {
	set wave2 [showWave $optic $i]
	set refract [diffRefract $wave2 $wave1 $zenith]
	setFocal $optic $i yoff [expr -1.*$refract/60.]
	}
echo Max refraction [format %.2f $refract]
   return
   }

################################################################
#Assume I have run prismSet.  Solve for optimal prism angles.
#I ONLY use the center point.
proc prismSolve {optic} {

#Clean out any old least-squares flags.
   clearFlags

#Clean out any old spots.
   clearSpot
   addSpot 0 0
#   addSpot 0 1
#   addSpot 0 -1
#   addSpot 1 0
#   addSpot -1 0
   set ncolor [exprGet $optic.ncolor]
   set colors ""

#Get list of colors
   for {set i 1} {$i <= $ncolor} {incr i} {
	lappend colors $i
	}

#Initial calculation of phi.
   set refract [diffRefract .55 .39 $zenith]
   set phi [expr ($refract/1.6)*.5]

#Use an on-axis spot
   set ncolor [exprGet $optic.ncolor]
   set nsurf [exprGet $optic.ncolor]

#Rotate prisms.  Run a sample ray and find prisms that way
   set surfids [surfGet $optic 1]
   set adc1 ""
   set adc2 ""
   foreach surfid $surfids {
	set name [string tolower [showName $optic $surfid]]
	if {$name == "adc1"} {lappend adc1 $surfid}
	if {$name == "adc2"} {lappend adc2 $surfid}
	}

   foreach surfid $adc1 {
	setSurf $optic $surfid phi [expr -1.*$phi]
	setSurfInc $optic $surfid phi .001
	}
   foreach surfid $adc2 {
	setSurf $optic $surfid phi $phi
	setSurfInc $optic $surfid phi .001
	}

#Set inc to a fixed amount
   if {[llength $adc1] == 0 || [llength $adc2] == 0} return
   linkColorFlag $optic $colors 2
   linkSurfFlag $optic $adc1 phi 4
   linkSurfFlag $optic $adc2 phi 6

   global _optic
   set _optic(lscale) 0
   set _optic(niter) 3
   set _optic(fractmax) 1

   set ang1 [format %.2f [showSurf $optic [lindex $adc1 0] phi]]
   set ang2 [format %.2f [showSurf $optic [lindex $adc2 0] phi]]

   pfastLstsq $optic

   set ang1 [format %.2f [showSurf $optic [lindex $adc1 0] phi]]
   set ang2 [format %.2f [showSurf $optic [lindex $adc2 0] phi]]
   echo phi1 = $ang1, phi2 = $ang2
   return
   }

##########################################################################
#Reset prisms 
proc prismReset {optic} {

#Get list of adc surfaces.
   set surfids [surfGet $optic 1]
   set adc1 ""
   set adc2 ""
   foreach surfid $surfids {
	set name [string tolower [showName $optic $surfid]]
	if {$name == "adc1"} {setSurf $optic $surfid phi 0}
	if {$name == "adc2"} {setSurf $optic $surfid phi 0}
	}
   set ncolor [exprGet $optic.ncolor]
   for {set i 1} {$i <= $ncolor} {incr i} {
	setFocal $optic $i yoff 0
	}
   return
   }

##########################################################################
#Evaluate adc for a range of angles for a specified set of filters.
#icolor is used as a template configuration.
#This routine works even for designs without adc (like may11).

proc adcEval {inoptic rad icolor series filter} {
   set optic [filterSample $inoptic $icolor $series $filter]
   set ncolor [exprGet $optic.ncolor]
   set colors ""
   for {set i 1} {$i <= $ncolor} {incr i} {
	lappend colors $i
	}
   foreach zenith "0 30 45 55 60" {
	prismSet $optic $zenith
	if {$zenith > 0} {prismSolve $optic}
	set avg 0
	foreach i "-1 1" {
	   set avg [expr $avg + [raySum $optic [expr $i*$rad] 0 $colors]]
	   }
	foreach i "-1 1" {
	   set avg [expr $avg + [raySum $optic 0 [expr $i*$rad] $colors]]
	   }
	set d80($zenith) [expr $avg/4.]
#	spotPlot $optic $xmm $ymm $colors
	prismReset $optic
	}
   echo
   foreach zenith "0 30 45 55 60" {
	echo Zenith angle $zenith, radius $rad, D80= [format %.2f \
	   $d80($zenith)]
 	}
   opticDel $optic
   return
   }

##########################################################################
#Simulate images at different zenith angles.

proc weakSim {inoptic xmm ymm icolor series filter zeniths} {

   set optic [filterSample $inoptic $icolor $series $filter]
   set ncolor [exprGet $optic.ncolor]
   set colors ""
   for {set i 1} {$i <= $ncolor} {incr i} {
	lappend colors $i
	}

#Come up with a useful name for the output file
   set name [lindex [lindex [exprGet $optic.name] 0] 0]
   set name [file root $name]

   global SIGMA NPIX PIXMM
   set NPIX 64
   set PIXMM .015
   foreach zenith $zeniths {
	set xname $xmm
	regsub -- {-} $xname M xname
	set yname $ymm
	regsub -- {-} $yname M yname
	set file $name-x${xname}-y${yname}-$filter-z$zenith.fit
	prismSet $optic $zenith
	if {$zenith > 0} {prismSolve $optic}

#SIGMA is .6 arcsec at zenith for .6 microns
	set SIGMA [expr .7/2.35]
#	set SIGMA [expr ($SIGMA/pow(cos($zenith/60.),.6))]
	set wave [showWave $optic 1]
#	set SIGMA [expr $SIGMA*pow(.6/$wave,.2)]
	spotMap $optic $xmm $ymm $colors
	exec mv spot.fit $file
	}	
   opticDel $optic
   return
   }

#########################################################################
#Top level job to process multiple designs, filters, etc.
#Things are a bit hard-wired for my filter scheme.

proc weakSimJob {hndls positions series filters zeniths} {
   set filtername(g) gprime
   set filtername(r) rprime
   set filtername(i) iprime
   set filtername(z) zprime
   set filternum(g) 1
   set filternum(r) 3
   set filternum(i) 5
   set filternum(z) 7
   foreach hndl $hndls {
	foreach filter $filters {
	   foreach position $positions {
		set xmm [lindex $position 0]
		set ymm [lindex $position 1]
		weakSim $hndl $xmm $ymm $filternum($filter) $series \
		   $filtername($filter) $zeniths
		}
	   }
	}
   return
   }

#########################################################################
#The following runs in kentools.
#In cray, I create an over-sampled image.  I should have one oversampled
#by 10x relative to the biggest desired pixel size.
#	set NPIX 256
#	set SIGMA [expr .6/2.35]
#	set PIXMM [expr .015/8.]
#	spotMap hndl xmm ymm colors
#Here, I read it, shift it to a set of different pixel centers,
#then compress it down to the expected Decam pixel scale.

#Allow for a gain for adding photon statistics
proc weakSample {file {gain 10.}} {
   set ext [file extension $file]
   if {$ext == ""} {
	set file $file.fit
	}
   set reg [regNew]

   regReadAsFits $reg $file

#Get pixel scale
   set scale [hdrGetAsDbl $reg.hdr SCALE]

#Shift and subsample does not appear to make any difference, so I won't
#try it here.  Instead, test difference compression factors
   set npix [exprGet $reg.nrow]
   set fwhmList ""
   set scaleList ""
   loop i 1 11 {
	set avg($i) 0.
	}
   set ncount 0
   loop count 0 20 {
   puts stdout "  i  Scale    FWHM  Whisk    PA  SymWhisk"
   loop i 1 11 {
	set new [regNewFromReg $reg]
	imgPhoton $new $gain

#Shift to new center.  I am oversampled by a factor 8, so use that as my
#shift window.
	set dx [expr round(8*([random 1000]/1000.-.5))]
	set dy [expr round(8*([random 1000]/1000.-.5))]
	regIshift $new $dx $dy
	set comp [imgCompress $new $i]
	set r [expr round($npix/(2.*$i))]
	axes $comp $r $r 128 128
	global img
	set hndl $img(axes)
	set sigmaj [exprGet $hndl.sigmaj]
	set sigmin [exprGet $hndl.sigmin]
	set len [expr sqrt($sigmaj*$sigmaj - $sigmin*$sigmin)]
	set fwhm [expr sqrt($sigmaj*$sigmin)]
	set pixscale [expr $scale*1.*$i]
	set len [expr $len*$pixscale*2.35]
	set fwhm [expr $fwhm*$pixscale*2.35]
	if {$i == 1} {set fwhm0 $fwhm}
	set symwhisk [expr sqrt(max(0.,$fwhm*$fwhm - $fwhm0*$fwhm0))]
	set ang [exprGet $hndl.pos_maj]
	puts stdout [format "%3d %6.3f %6.3f %6.3f %6.3f %6.3f" $i \
		$pixscale $fwhm $len $ang $symwhisk]
	lappend fwhmList $fwhm
	lappend scaleList $pixscale
	regDel $new
	regDel $comp
	set avg($i) [expr $avg($i) + $len]
	}
	incr ncount
	}
   puts stdout ""
   puts stdout "  i  Scale    Whisk"
   loop i 1 11 {
	set pixscale [expr $scale*1.*$i]
	set avg($i) [expr $avg($i)/$ncount]
	puts stdout [format "%3d %6.3f %6.3f" $i $pixscale $avg($i)]
	}
   regDel $reg
   pgEnv 0 .4 .6 .75 0 0
   pgLine $scaleList $fwhmList
   return
   }

#########################################################################
#The following routine runs in kentools.  It determines all fits files,
#reads in the images and computes axis info, then prints out a nice table.

proc weakAxes {} {
   set files [lsort [glob -nocomplain *.fit]]
   set n 0
   global pointer sigavg eps
   foreach var "pointer sigavg eps" {
	if {[info exists $var]} {unset $var}
	}

   foreach file $files {
	set reg [openimage $file]

#I've standardized on 64 pix file size.
	if {[exprGet $reg.nrow] != 64} {
	   regDel $reg
	   continue
	   }
	set list [split [file root $file] -]
	set design($n) [lindex $list 0]
	set x($n) [lindex $list 1]
	set y($n) [lindex $list 2]
	set filter($n) [lindex $list 3]
	set zenith($n) [lindex $list 4]
	set index $design($n),$x($n),$y($n),$filter($n),$zenith($n)
	set pointer($index) $n
	axes $reg 32 32
	global img
	set hndl $img(axes)
	set eps($n) [format %.4f [exprGet $hndl.eps]]
	set sigmaj [exprGet $hndl.sigmaj]
	set sigmin [exprGet $hndl.sigmin]
	set sigavg($n) [format %.2f [expr ($sigmaj+$sigmin)/2.]]
	incr n
	regDel $reg
	}
   set format "%-6s %-4s %-4s %-6s %-4s %.2f   %.4f"
   puts stdout "Design x    y    filter zang sigavg   eps"
   loop i 0 $n {
	puts stdout [format $format $design($i) $x($i) $y($i) $filter($i) \
	   $zenith($i) $sigavg($i) $eps($i)]
	}
   return
   }

#####################################################################
proc weakList {} {
   global pointer sigavg eps
   set format1 "%5s %3s    %8s       %8s       %8s       %8s     "
   set format1a "%5s %3s  %6s%7s  %6s%7s  %6s%7s  %6s%7s"
   set format2 "%-5s %3s %6.2f%7.3f  %6.2f%7.3f  %6.2f%7.3f  %6.2f%7.3f"

   puts stdout [format $format1 Filt zen may11 may11a adc5 vista]
   puts stdout [format $format1a "" "" FWHM eps FWHM eps FWHM eps FWHM eps]
   foreach filter "gprime rprime iprime zprime" {
	foreach zenith "z0 z35 z55" {
	   set scale(may11) 17.65
	   set scale(may11a) 17.65
	   set scale(adc5) 17.71
	   set scale(vista) 16.91
	   foreach design "may11 may11a adc5 vista" {
		set sigsum($design,$filter,$zenith) 0.
		set epssum($design,$filter,$zenith) 0.
		foreach pos "xM225,y0 x225,y0 x0,y225 xM225,y0" {
		   set i $pointer($design,$pos,$filter,$zenith)
		   set sigsum($design,$filter,$zenith) [expr \
			$sigsum($design,$filter,$zenith) + $sigavg($i)]
		   set epssum($design,$filter,$zenith) [expr \
			$epssum($design,$filter,$zenith) + $eps($i)]
		   }

#Convert to FWHM in arcsec.
		set sigsum($design,$filter,$zenith) [expr $scale($design) * \
	.015 * $sigsum($design,$filter,$zenith) *2.35/ 4.]
		set epssum($design,$filter,$zenith) [expr \
	$epssum($design,$filter,$zenith) / 4.]
		}
	   }
	}
   foreach filter "gprime rprime iprime zprime" {
	regsub prime $filter ' filname
	foreach zenith "z0 z35 z55" {
	   regsub z $zenith "" zenname
	   puts stdout [format $format2 $filname $zenname \
		$sigsum(may11,$filter,$zenith) $epssum(may11,$filter,$zenith) \
		$sigsum(may11a,$filter,$zenith) $epssum(may11a,$filter,$zenith) \
		$sigsum(adc5,$filter,$zenith) $epssum(adc5,$filter,$zenith) \
		$sigsum(vista,$filter,$zenith) $epssum(vista,$filter,$zenith)]
	   }
	puts stdout ""
	}
   return
   }
#####################################################################
#List the broadening function = eps*sig,  not just sig or eps.
proc weak2List {} {
   global pointer sigavg eps
   set format1 "%-5s %3s %8s %8s %8s %8s"
   set format2 "%-5s %3s %8.3f %8.3f %8.3f %8.3f"

   puts stdout [format $format1 Filt zen may11 may11a adc5 vista]
   foreach filter "gprime rprime iprime zprime" {
	foreach zenith "z0 z35 z55" {
	   set scale(may11) 17.65
	   set scale(may11a) 17.65
	   set scale(adc5) 17.71
	   set scale(vista) 16.91
	   foreach design "may11 may11a adc5 vista" {
		set sigsum($design,$filter,$zenith) 0.
		set epssum($design,$filter,$zenith) 0.
		foreach pos "xM225,y0 x225,y0 x0,y225 xM225,y0" {
		   set i $pointer($design,$pos,$filter,$zenith)
		   set sigsum($design,$filter,$zenith) [expr \
			$sigsum($design,$filter,$zenith) + $sigavg($i)]
		   set epssum($design,$filter,$zenith) [expr \
			$epssum($design,$filter,$zenith) + $eps($i)]
		   }

#Convert to FWHM in arcsec.
		set sigsum($design,$filter,$zenith) [expr $scale($design) * \
	.015 * $sigsum($design,$filter,$zenith) *2.35/ 4.]
		set epssum($design,$filter,$zenith) [expr \
	$epssum($design,$filter,$zenith) / 4.]
		}
	   }
	}
   foreach filter "gprime rprime iprime zprime" {
	regsub prime $filter ' filname
	foreach zenith "z0 z35 z55" {
	   regsub z $zenith "" zenname
	   set broad1 [expr $sigsum(may11,$filter,$zenith) * \
		$epssum(may11,$filter,$zenith)]
	   set broad2 [expr $sigsum(may11a,$filter,$zenith)  * \
		$epssum(may11a,$filter,$zenith)]
	   set broad3 [expr $sigsum(adc5,$filter,$zenith) * \
		$epssum(adc5,$filter,$zenith)]
	   set broad4 [expr $sigsum(vista,$filter,$zenith) * \
		$epssum(vista,$filter,$zenith)]
	   puts stdout [format $format2 $filname $zenname \
		$broad1 $broad2 $broad3 $broad4]
	   }
	puts stdout ""
	}
   return
   }

#####################################################################
#Amplitude difference between two "whiskers"
#x0, y0 are vector components of first whisker.
#x1, y1 etc. for second.

proc whiskDiff {x0 y0 x1 y1} {

#Difference in amplitude is tricky - go back to matrix of 2nd moments
   set dxx [expr $x1*$x1 - $x0*$x0]
   set dyy [expr $y1*$y1 - $y0*$y0] 
   set dxy [expr $x1*$y1 - $x0*$y0]
   set ds [expr sqrt(sqrt(pow($dxx-$dyy,2) + pow(2.*$dxy,2)))]
   return $ds
   }

#####################################################################
#Amplitude difference between two polarization vectors.
#e11, e12 are vector components of first whisker.
#e21, e22 etc. for second.

proc polDiff {e11 e12 e21 e22} {

#Difference in amplitude is tricky - go back to matrix of 2nd moments
   set dxx [expr $x1*$x1 - $x0*$x0]
   set dyy [expr $y1*$y1 - $y0*$y0] 
   set dxy [expr $x1*$y1 - $x0*$y0]
   set ds [expr sqrt(sqrt(pow($e21-$e11,2) + pow($e22 - $e12,2)))]
   return $ds
   }

#####################################################################
#Amplitude difference between two "dilutions"

proc diluteDiff {amp1 amp2} {

#Difference in amplitude is tricky - go back to matrix of 2nd moments
   set ds [expr sqrt(sqrt(pow( pow($amp1,2) - pow($amp2,2) ,2) ))]
   return $ds
   }

#####################################################################
#Amplitude and position angle difference between two "whiskers"
#x0, y0 are vector components of first whisker.
#x1, y1 etc. for second.
#If I swap x0,y0 and x1,y1, the positon angle change by 90 deg.
#I can use this proc to make "whisker" plots.

proc whiskDiffPA {x0 y0 x1 y1} {

#Difference in amplitude is tricky - go back to matrix of 2nd moments
   set dxx [expr $x1*$x1 - $x0*$x0]
   set dyy [expr $y1*$y1 - $y0*$y0] 
   set dxy [expr $x1*$y1 - $x0*$y0]
   set ds [expr sqrt(sqrt(pow($dxx-$dyy,2) + pow(2.*$dxy,2)))]

#Position angle.  Taken from "whisker"
   set s2t [expr 2.*$dxy]
   set c2t [expr $dxx-$dyy]
   set theta [expr 180.*atan2($s2t,$c2t)/(2.*3.141593)]
   return [list $ds $theta]
   }

#####################################################################
#Compute quantities for benefit of PPARC proposal.  Just do one design.

proc pparc {hndl} {
   set seeing .7
   global pointer sigavg eps
   set filtername(g) gprime
   set filtername(r) rprime
   set filtername(i) iprime
   set filtername(z) zprime
   set filternum(g) 1
   set filternum(r) 3
   set filternum(i) 5
   set filternum(z) 7
   global fwhm fwhm0 xwhm ywhm ang

#Feed in sky angles in degrees
   set thetalist "0 0.5 1.0"
   set nang [llength $thetalist]
   set poslist [list [list 1 0] [list -1 0] [list 0 1] [list 0 -1]]
   set npos [llength $poslist]
   set zeniths "0 30 35 45 55 60"
   set filters "r i"
   if {![info exists fwhm]} {
	set design [showDesign $hndl]
	foreach filter $filters {
	   set filterName $filtername($filter)
	   set filterNum $filternum($filter)
	   foreach zenith $zeniths {
		set optic [filterSample $hndl $filterNum vista $filterName]
		prismSet $optic $zenith
		if {$zenith > 0} {prismSolve $optic}
		set ncolor [exprGet $optic.ncolor]
		set colors ""
		for {set i 1} {$i <= $ncolor} {incr i} {lappend colors $i}

#Set wavelengths and ADC position.
		set iang 0
		foreach theta $thetalist {
		   set theta [expr $theta*60.]
		   set list [skytofocal $optic $theta 0 1]
		   set rad [lindex $list 0]
		   set ipos 0
		   foreach pos $poslist {
			set xmm [expr [lindex $pos 0] * $rad]
			set ymm [expr [lindex $pos 1] * $rad]
			set d80 [raySum $optic $xmm $ymm $colors]
#Approx. dilution derived from d80
			set fwhm0($design,$filter,$zenith,$iang,$ipos) \
			   [expr $d80/1.53]
			set list [whisker $optic $xmm $ymm $colors]
			set fwhm($design,$filter,$zenith,$iang,$ipos) \
			   [lindex $list 2]
			set ang($design,$filter,$zenith,$iang,$ipos) \
			   [lindex $list 3]
			incr ipos
			}
		   incr iang
		   }
		opticDel $optic
		}
	   }
	}

#At each zenith angle, compute ellipticity for a specific seeing.
#vector length of the whisker
   set design [showDesign $hndl]
   foreach filter $filters {
	set filterName $filtername($filter)
	set filterNum $filternum($filter)
	foreach zenith $zeniths {
	   loop iang 0 $nang {
		set epssum 0.
		set fwhmsum 0.
		loop ipos 0 $npos {
		   set w0 $fwhm($design,$filter,$zenith,$iang,$ipos)
		   set w $fwhm($design,$filter,$zenith,$iang,$ipos)
		   set b [expr sqrt($seeing*$seeing + $w0*$w0 - $w*$w/2.)]
		   set a [expr sqrt($b*$b + $w*$w)]
		   set eps [expr ($a-$b)/$a]
		   set fwhmsum [expr $fwhmsum + sqrt(($a*$a + $b*$b)/2.)]
		   set epssum [expr $epssum + $eps]
		   }
		set epsavg($design,$filter,$zenith,$iang) \
		   [expr $epssum/$npos]
		set fwhmavg($design,$filter,$zenith,$iang) \
		   [expr $fwhmsum/$npos]
		}
	   }
	}


#Print out results.
   set format1 "%8s %8s %8s %8s %8s %8s"
   set format2 "%8s %8s %8.0f %8.1f %8.3f %8.3f"
   foreach filter $filters {
	puts stdout [format $format1 design filter zenith angle fwhm eps]
	foreach zenith $zeniths {
	   loop iang 0 $nang {
		set theta [lindex $thetalist $iang]
		puts stdout [format $format2 $design $filter $zenith $theta \
		   $fwhmavg($design,$filter,$zenith,$iang) \
		   $epsavg($design,$filter,$zenith,$iang)]
		}
	   puts stdout ""
	   }
	puts stdout ""
	}
   return
   }

##################################################################
#Computes change in PSF shape over focal plane due to defocus
#I will input handles of the same design with multiple focii.
#Note:  FWHM is the whisker length, not the dilution!!!  I hae been telling
#people that it is the dilution, incorrectly.

proc whiskFocus {hndls} {
   set filtername(g) gprime
   set filtername(r) rprime
   set filtername(i) iprime
   set filtername(z) zprime
   set filternum(g) 1
   set filternum(r) 3
   set filternum(i) 5
   set filternum(z) 7
   set format1 "%5s %3s    %8s"
   set format1a "%5s %3s  %6s%7s"
   set format2 "%-5s %3s %6.2f%7.2f"

   global ffwhm fxwhm fywhm fang

   set poslist [list [list -225 0] [list -160 0] [list 0 0] [list 160 0] \
	[list 225 0] [list 0 -225] [list 0 -160] [list 0 160] [list 0 225]]
   set npos [llength $poslist]

if {![info exists ffwhm]} {
   loop k 0 [llength $hndls] {
	set hndl [lindex $hndls $k]
	set design $k
	foreach filter "g r i z" {
	   set filterName $filtername($filter)
	   set filterNum $filternum($filter)
	   foreach zenith "0" {
		set optic [filterSample $hndl $filterNum vista $filterName]
		prismSet $optic $zenith
		if {$zenith > 0} {prismSolve $optic}
		set ncolor [exprGet $optic.ncolor]
		set colors ""
		for {set i 1} {$i <= $ncolor} {incr i} {lappend colors $i}

#Set wavelengths and ADC position.
		set ipos 0
		foreach pos $poslist {
		   set xmm [lindex $pos 0]
		   set ymm [lindex $pos 1]
		   set list [whisker $optic $xmm $ymm $colors]
		   set fxwhm($design,$filter,$zenith,$ipos) [lindex $list 0]
		   set fywhm($design,$filter,$zenith,$ipos) [lindex $list 1]
		   set ffwhm($design,$filter,$zenith,$ipos) [lindex $list 2]
		   set fang($design,$filter,$zenith,$ipos) [lindex $list 3]
		   incr ipos
		   }
		opticDel $optic
		}
	   }
	}
}

#At each zenith angle, compute the average FWHM and the average change in
#vector length of the whisker from differencing positions.
   set name [showDesign [lindex $hndls 0]]
   puts stdout [format $format1 Filt zen $name]
   puts stdout [format $format1a "" "" FWHM diff]
   foreach filter "g r i z" {
	set filterName $filtername($filter)
	set filterNum $filternum($filter)
	foreach zenith "0" {
	   set a 0.
	   set d 0.
	   set n 0
	   loop ipos 0 $npos {
		loop k 1 [llength $hndls] {
		   set k0 [expr $k-1]

#Difference in amplitude is tricky - go back to matrix of 2nd moments
		   set x1 $fxwhm($k,$filter,$zenith,$ipos)
		   set x0 $fxwhm($k0,$filter,$zenith,$ipos)
		   set y1 $fywhm($k,$filter,$zenith,$ipos)
		   set y0 $fywhm($k0,$filter,$zenith,$ipos)
		   set ds [whiskDiff $x0 $y0 $x1 $y1]
		   set d [expr $d + $ds]
		   set a [expr $a + \
			($ffwhm($k,$filter,$zenith,$ipos) + \
			$ffwhm($k0,$filter,$zenith,$ipos))/2.]
		   incr n
		   }
		}

#The average is slightly wrong.
	   set fwhmavg($filter,$zenith) [expr $a/$n]
	   set fwhmdiff($filter,$zenith) [expr $d/$n]
	   }
	}

#Print out results.
   foreach filter "g r i z" {
	foreach zenith "0" {
	   puts stdout [format $format2 $filter $zenith \
		$fwhmavg($filter,$zenith) \
		$fwhmdiff($filter,$zenith)]
	   }
	puts stdout ""
	}
   return
   }

##################################################################
#Plot the PSF shape over the focal plane.
#Specify the filter, because I don't necessarily case about all filters

proc whiskMap {hndl ifil} {
   plotInit a
   set xrad [showFocal $hndl $ifil xrad]
   set yrad [showFocal $hndl $ifil yrad]
   set scale [expr abs((.01*2.*$xrad/.1))]

#For weirdtrend image, make bigger lines
#   set scale [expr $scale*2.]
   if {$xrad < 0} {
	set xmin $xrad
	set xmax [expr -1.*$xrad]
   } else {
	set xmin -$xrad
	set xmax $xrad
	}
   if {$yrad < 0} {
	set ymin $yrad
	set ymax [expr -1.*$yrad]
   } else {
	set ymin -$yrad
	set ymax $yrad
	}
   pgEnv $xmin $xmax $ymin $ymax 1 0
   pgLabel X Y "PSF Whisker Map"
   pgSfs 2
   set ampmax 0.
   set dilutemax 0.
   set amprms 0.
   set nrms 0
   set nstep 10
   loop i 0 $nstep {
	set xmm [expr $xrad*($i-($nstep/2. - 0.5))/($nstep/2.)]
	loop j 0 $nstep {
	   set ymm [expr $yrad*($j-($nstep/2. - 0.5))/($nstep/2.)]
	   if {$xrad > 0 && $yrad > 0} {
		if {pow($xmm/$xrad,2) + pow($ymm/$yrad,2) > 1.} continue
		}
	   set list [whisker $hndl $xmm $ymm $ifil]
	   set amp [lindex $list 2]
	   set pa [lindex $list 3]
	   set ampmax [expr max($amp,$ampmax)]
	   set amprms [expr $amprms + $amp*$amp]
	   incr nrms
	   set amp [expr $amp*$scale]
	   set ang [expr $pa/57.3]
	   set x1 [expr $xmm - $amp*cos($ang)/2.]
	   set x2 [expr $xmm + $amp*cos($ang)/2.]
	   set y1 [expr $ymm - $amp*sin($ang)/2.]
	   set y2 [expr $ymm + $amp*sin($ang)/2.]
	   pgSci 1
	   pgLine "$x1 $x2" "$y1 $y2"

	   pgSci 1
	   update
	   }
	}
#Draw a line 1"
   set xmin [expr -.5*$scale]
   set xmax [expr .5*$scale]
   pgSci 7
   pgLine "$xmin $xmax" "0 0"
   pgSci 1
   echo [format "Max. whisker length %.2f" $ampmax] arcsec
   set amprms [expr sqrt($amprms/$nrms)]
   echo [format "rms  whisker length %.2f" $amprms] arcsec
   return
   }

##################################################################
#Plot change in PSF shape between 2 designs.
#The 2 designs should have all parameter differences in them already.
#Specify the filter, because I don't necessarily case about all filters

proc whiskPlot {hndl1 hndl2 ifil} {
   plotInit a
   set xrad [showFocal $hndl1 $ifil xrad]
   set yrad [showFocal $hndl1 $ifil yrad]
   set scale [expr (.01*2.*$xrad/.1)]
   pgEnv -$xrad $xrad -$yrad $yrad 1 0
   pgLabel X Y "PSF Dilution and Polarization changes"
   pgSfs 2
   set ampmax 0
   set dilutemax 0
   set e1avg 0.
   set e2avg 0.
   set n 0
   loop i 0 10 {
	set xmm [expr $xrad*($i-4.5)/5.]
	loop j 0 10 {
	   set ymm [expr $yrad*($j-4.5)/5.]
	   if {$xrad > 0 && $yrad > 0} {
		if {pow($xmm/$xrad,2) + pow($ymm/$yrad,2) > 1.} continue
		}
	   set list [polarize $hndl1 $xmm $ymm $ifil]
	   set e11 [lindex $list 0]
	   set e12 [lindex $list 1]
	   set e13 [lindex $list 2]
	   set list [polarize $hndl2 $xmm $ymm $ifil]
	   set e21 [lindex $list 0]
	   set e22 [lindex $list 1]
	   set e23 [lindex $list 2]
	   set e1($i,$j) [expr $e21-$e11]
	   set e2($i,$j) [expr $e22-$e12]
	   set e3($i,$j) [expr abs($e23-$e13)]

	   set e1avg [expr $e1avg + $e1($i,$j)]
	   set e2avg [expr $e2avg + $e2($i,$j)]
	   incr n
	   }
	}
   set e1avg [expr $e1avg/$n]
   set e2avg [expr $e2avg/$n]
   set ampavg 0.
   loop i 0 10 {
	set xmm [expr $xrad*($i-4.5)/5.]
	loop j 0 10 {
	   set ymm [expr $yrad*($j-4.5)/5.]
	   if {$xrad > 0 && $yrad > 0} {
		if {pow($xmm/$xrad,2) + pow($ymm/$yrad,2) > 1.} continue
		}
	   set list [polarize $hndl1 $xmm $ymm $ifil]
	   set e11 [lindex $list 0]
	   set e12 [lindex $list 1]
	   set e13 [lindex $list 2]
	   set list [polarize $hndl2 $xmm $ymm $ifil]
	   set e21 [lindex $list 0]
	   set e22 [lindex $list 1]
	   set e23 [lindex $list 2]
	   set amp [expr sqrt(sqrt(pow($e1($i,$j) - $e1avg,2) + \
		pow($e2($i,$j) - $e2avg,2)))]
	   set ang [expr atan2($e2($i,$j) - $e2avg, $e1($i,$j) - $e1avg)/2.]
	   set ampmax [expr max($amp,$ampmax)]
	   set ampavg [expr $ampavg + $amp]

#Scale for plotting.
	   set amp [expr $amp*$scale]
	   set x1 [expr $xmm - $amp*cos($ang)/2.]
	   set x2 [expr $xmm + $amp*cos($ang)/2.]
	   set y1 [expr $ymm - $amp*sin($ang)/2.]
	   set y2 [expr $ymm + $amp*sin($ang)/2.]
	   pgSci 1
	   pgLine "$x1 $x2" "$y1 $y2"

#Also get dilution difference
	   set diff [expr sqrt($e3($i,$j))]
	   set dilutemax [expr max($diff,$dilutemax)]
	   pgSci 1
	   update
	   }
	}

#Draw a line 1"
   set xmin [expr -.5*$scale]
   set xmax [expr .5*$scale]
   pgSci 7
   pgLine "$xmin $xmax" "0 0"
   pgSci 1
   set ampavg [expr $ampavg/$n]
   echo [format "Max. whisker length %.2f" $ampmax] arcsec
   echo [format "AVg. whisker length %.2f" $ampavg] arcsec
   echo [format "Max. dilution length %.2f" $dilutemax] arcsec
   return
   }

##################################################################
#Make a "library" of whisker plots.
#I will read a set of adjustments to make.
#I will make a whisker different plot between input design and the tweaked
#design, the tweaked design and the tweaked design with a defocus,
#and the input design and the tweaked design with defocus.
#Tweaks are read from a file "tweak"

proc whiskLib {hndl0 ifil} {
   if {![file exists tweak]} {
	echo Need a file \"tweek\"!
	return
	}
   set fid [open tweak]
   while {1} {
	set line [string trim [gets $fid]]
	if {[eof $fid]} break
	if {[string length $line] == 0} continue
	if {[string index $line 0] == "#"} continue

#Line is a descriptive name followed by a TCL instruction with % in place
#of the optic structure handle.
	set name [lvarpop line]
	set hndl1 [opticNew]
	set hndl2 [opticNew]
	set hndl3 [opticNew]
	opticCopy $hndl0 $hndl1
	opticCopy $hndl0 $hndl2
	opticCopy $hndl0 $hndl3
	regsub {%} $line $hndl1 line1
	regsub {%} $line $hndl2 line2
	regsub {%} $line $hndl3 line3
	eval $line1
	eval $line2

#Hardwired: Focus is surface 2.
	incrSurf $hndl2 2 z -.1
	incrSurf $hndl3 2 z .1

	set gif1 $name-1.gif
	set gif2 $name-2.gif
	set gif3 $name-3.gif
	if {![file exists $gif1]} {
	   pgBegin $gif1/GIF
	   pgGeomSet 400 400 0 0
	   whiskPlot $hndl0 $hndl1 $ifil
	   m2
	   pgEnd
	   }

	if {![file exists $gif2]} {
	   pgBegin $gif2/GIF
	   pgGeomSet 400 400 0 0
	   whiskPlot $hndl0 $hndl2 $ifil
	   m2
	   pgEnd
	   }

	if {![file exists $gif3]} {
	   pgBegin $gif3/GIF
	   pgGeomSet 400 400 0 0
	   whiskPlot $hndl0 $hndl3 $ifil
	   m2
	   pgEnd
	   }

	lappend names $name
	opticDel $hndl1
	opticDel $hndl2
	opticDel $hndl3
	}
   close $fid

#Create html file
   set design [showDesign $hndl0]
   set fid [open tweak.html w]
   puts $fid "<head>"
   puts $fid "<title>Design tweaks</title>"
   puts $fid "</head>"
   puts $fid "<body bgcolor=white>"
   puts $fid "<h1 align=center>	Design tweaks for $design</h1>"
   puts $fid "<table align=center>"
   puts $fid "<tr><th>Name</th><th>Tweak</th>"
   puts $fid "<th>Tweak plus -Defocus</th>"
   puts $fid "<th>Tweak plus +Defocus</th></tr>"
   foreach name $names {
	puts $fid "<tr align=center><td>$name</td>"
	puts $fid "<td><img src=\"$name-1.gif\"></td>"
	puts $fid "<td><img src=\"$name-2.gif\"></td>"
	puts $fid "<td><img src=\"$name-3.gif\"></td></tr>"
	puts $fid ""
	}
   puts $fid "</table></body>"
   close $fid
   return
   }

##################################################################
#Computes change in PSF shape over focal plane due to defocus
#I will input handles of the same design with multiple focii.
#In this version, I generalize things.

proc lsstFocus {hndls filters} {
   set format1 "%5s %3s    %8s"
   set format1a "%5s %3s  %6s%7s"
   set format2 "%-5s %3s %6.2f%7.2f"

   set hndl0 [lindex $hndls 0]
   set filter0 [lindex $filters 0]
   set rad [expr abs([showFocal $hndl0 $filter0 xrad])]
   set poslist ""
   foreach x "-1 -.7 0 .7 1" {
	lappend poslist "[expr $x*$rad] 0"
	}
   foreach y "-1 -.7 .7 1" {
	lappend poslist "0 [expr $y*$rad]"
	}
   set npos [llength $poslist]

   set zenith 0
   loop k 0 [llength $hndls] {
	set optic [lindex $hndls $k]
	set design $k

#Set wavelengths and ADC position.
	set ipos 0
	foreach pos $poslist {
	   set xmm [lindex $pos 0]
	   set ymm [lindex $pos 1]
	   set list [whisker $optic $xmm $ymm $filters]
	   set fxwhm($design,$filters,$zenith,$ipos) [lindex $list 0]
	   set fywhm($design,$filters,$zenith,$ipos) [lindex $list 1]
	   set ffwhm($design,$filters,$zenith,$ipos) [lindex $list 2]
	   set fang($design,$filters,$zenith,$ipos) [lindex $list 3]
	   incr ipos
	   }
	}

#At each zenith angle, compute the average FWHM and the average change in
#vector length of the whisker from differencing positions.
   set name [showDesign [lindex $hndls 0]]
   puts stdout [format $format1 Filt zen $name]
   puts stdout [format $format1a "" "" FWHM diff]
   set a 0.
   set d 0.
   set n 0
   loop ipos 0 $npos {
	loop k 1 [llength $hndls] {
	   set k0 [expr $k-1]

#Difference in amplitude is tricky - go back to matrix of 2nd moments
	   set x1 $fxwhm($k,$filters,$zenith,$ipos)
	   set x0 $fxwhm($k0,$filters,$zenith,$ipos)
	   set y1 $fywhm($k,$filters,$zenith,$ipos)
	   set y0 $fywhm($k0,$filters,$zenith,$ipos)
	   set ds [whiskDiff $x0 $y0 $x1 $y1]
	   set d [expr $d + $ds]

#The average is wrong because I leave out handle 0.  OK, now average it in.
	   set a [expr $a + \
		($ffwhm($k,$filters,$zenith,$ipos) + \
		$ffwhm($k0,$filters,$zenith,$ipos))/2.]
	   incr n
	   }
	}
   set fwhmavg($filters,$zenith) [expr $a/$n]
   set fwhmdiff($filters,$zenith) [expr $d/$n]

#Print out results.
   puts stdout [format $format2 $filters $zenith \
	$fwhmavg($filters,$zenith) \
	$fwhmdiff($filters,$zenith)]
   return
   }

#####################################################################
#Evaluate PF coma of blanco

proc pfComa {{sign 0}} {
   set fid [open pfcoma.txt]
   plotInit
   pgSci 1
#   pgEnv 0 75 0 1500 0 0
   pgEnv -1200 1200 -1200 1200 0 0
#   pgLabel "zenith ang" "Coma" "Blanco decenter"
   pgLabel X Y "Blanco decenter"
   set icol 1
   while {1} {
	incr icol
	if {$icol > 7} {
	   set icol 1
	   }
	pgSci $icol
	set line [gets $fid]
	if {[eof $fid]} break
	set az [lindex $line 0]
	loop i 0 5 {
	   set z($i) [lindex "0 15 30 45 60" $i]
	   set list [azelToEq $az [expr 90.-$z($i)] -32.]

#Position angle of zenith w.r.t. north
	   set pa($i) [keylget list pa]
	   set coma($i) [lindex $line [expr $i+1]]
	   set ang($i) [lindex $line [expr $i+6]]
	   set r [expr $ang($i)/57.3]
	   set x($i) [expr $coma($i)*cos($r)]
	   set y($i) [expr $coma($i)*sin($r)]
	   set xpa($i) [expr $z($i)*cos($pa($i)/57.3)]
	   set ypa($i) [expr $z($i)*sin($pa($i)/57.3)]
	   }
	echo
	echo Azimuth $az
	set xlist 0
	set ylist 0
	set xpalist 0
	set ypalist 0
	loop i 1 5 {
	   set x($i) [expr $x($i) - $x(0)]
	   set y($i) [expr $y($i) - $y(0)]
	   set mag [expr sqrt(pow($x($i),2) + pow($y($i),2))]
	   set phi [expr 57.3*atan2($y($i),$x($i))]
	   echo "   " z $z($i) mag [format %.0f $mag] PA [format %.0f $phi]

#Recompute x and y after adjusting phi from celestial to terrestrial
	   set phi [expr $phi + $sign*$pa($i)]
	   set x($i) [expr $mag*cos($phi/57.3)]
	   set y($i) [expr $mag*sin($phi/57.3)]
	   lappend xlist $x($i)
	   lappend ylist $y($i)

	   set xpa($i) [expr $xpa($i)*20.]
	   set ypa($i) [expr $ypa($i)*20.]
	   lappend xpalist $xpa($i)
	   lappend ypalist $ypa($i)
	   }
	pgSls 1
	pgLine $xlist $ylist
	pgSls 2
#	pgLine $xpalist $ypalist
	pgSls 1
#	pgText 62 $mag $az
	pgText $x(4) $y(4) $az
	}
   close $fid
   return
   }

#######################################################################
#Test dependence of aberration on offset from design focus
proc aberrCheck {optic focus focal ifil} {
   global _optic
   set _optic(fractmax) 1
   set _optic(niter) 1
   set hndl [opticNew]
   opticCopy $optic $hndl
   clearFlags
   setColorFlag $hndl $ifil 1
   setSurfFlag $hndl $focus z 1
   clearSpot
   addSpot 0 0
   addSpot 0 1
   set _optic(flagcache) $_optic(flags)
   set _optic(lscale) 0
   fastLstsq $hndl
   set edge [expr .8*[showFocal $hndl $ifil xrad]]

echo Nominal position
   zernikeMap3d $hndl 0 $edge $ifil

   incrSurf $hndl $focal z 2
   fastLstsq $hndl

echo Moved focal plane +2 mm
   zernikeMap3d $hndl 0 $edge $ifil

   incrSurf $hndl $focal z -4
   fastLstsq $hndl

echo Moved focal plane -2 mm
   zernikeMap3d $hndl 0 $edge $ifil
   opticDel $hndl
   return
   }

###################################################################
#Compare weak lensing folks definition with Hubble's.

proc eps {} {
   loop i 0 21 {
	set eps [expr $i/100.]
	set b/a [expr 1. - $eps]
	set ecc2 [expr 1.-pow([set b/a],2)]
	set weak [expr $ecc2/(2.-$ecc2)]
	echo Hubble [format %.2f $eps] Weak [format %.2f $weak]
	}
   return
   }

################################################################
#Check FWHM as we go through focus.  Assume we are at best focus.
#Do a set of 9 "exposures".  Look at center, edge
#Move the focal plane

proc focusSequence {hndl0 xmm ymm filter step} {
   set hndl [opticNew]

#Find focal plane
   set surfids [surfGet $hndl0 $filter]
   set focid -1
   foreach id $surfids {
	if {[regexp "primary" [string tolower [showName $hndl0 $id]]]} {
	   set focid $id
	   break
	   }
	}
   if {$focid < 0} {
	echo "No focal plane found!"
	return
	}
   echo Focal surface name $focid
   set fwhm0 [dilution $hndl0 $xmm $ymm $filter]

#Compute offset for minimum fwhm, ellipticity
   set fwhmmin 1000.
   set whiskmin 1000.
   foreach i "-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10" {
	opticCopy $hndl0 $hndl
	setSurf $hndl $focid z [expr [showSurf $hndl0 $focid z] + 1.*$i*$step]

#Center of field
	set fwhm [dilution $hndl $xmm $ymm $filter]

#Whisker
	set list [whisker $hndl $xmm $ymm $filter]
	set whisker [lindex $list 2]
	echo Frame $i FWHM (arcsec) [format %.3f $fwhm] Whisker (arcsec) \
	   [format %.3f $whisker]
	if {$fwhm < $fwhmmin} {
	   set ifwhmmin $i
	   set fwhmmin $fwhm
	   }
	if {$whisker < $whiskmin} {
	   set iwhiskmin $i
	   set whiskmin $whisker
	   }
	}
   opticDel $hndl
   return [list [format %.3f [expr 1.*$ifwhmmin*$step]] \
	[format %.3f [expr 1.*$iwhiskmin*$step]]]
   }

#########################################################################
#Determine difference between best focus and roundest images
#Step is in mm.

proc focSurf {hndl ifil step} {
   set xrad [showFocal $hndl $ifil xrad]
   set yrad [showFocal $hndl $ifil yrad]
   set circ 1
   if {$xrad < 0 && $yrad < 0} {set circ 0}

#For circular fields, just map the inscribed square.
   if {$circ} {
	set xrad [expr $xrad/sqrt(2.)]
	set yrad [expr $yrad/sqrt(2.)]
	}
   set diff() ""
   set xmin [expr -1.05*$xrad]
   set xmax [expr 1.05*$xrad]
   set ymin [expr -1.05*$yrad]
   set ymax [expr 1.05*$yrad]
   pgEnv $xmin $xmax $ymin $ymax 0 0
   foreach i "-2 -1 0 1 2" {
	foreach j "-2 -1 0 1 2" {
	   set xmm [format %.0f [expr $xrad*1.*$i/2.]]
	   set ymm [format %.0f [expr $yrad*1.*$j/2.]]
	   set list [focusSequence $hndl $xmm $ymm $ifil $step]
	   set fwhm [lindex $list 0]
	   set whisk [lindex $list 1]
	   set diff($xmm,$ymm) [format %.3f [expr $whisk-$fwhm]]
	   set diff($xmm,$ymm,xmm) $xmm
	   set diff($xmm,$ymm,ymm) $ymm
	   lappend diff() $xmm,$ymm
	   pgSch [expr abs($diff($xmm,$ymm))/.005]
	   if {$diff($xmm,$ymm) > 0} {
		pgPoint $xmm $ymm 2
	   } else {
		pgPoint $xmm $ymm 5
		}
	   }
	}
   foreach l $diff() {
	echo $diff($l,xmm) $diff($l,ymm) $diff($l)
	}
   return
   }

######################################################################
#Make a plot of the mosaic II focal plane.

proc m2 {} {
   pgSci 6
   pgRect -60 -30 -60 0
   pgRect -30 0 -60 0
   pgRect 0 30 -60 0
   pgRect 30 60 -60 0
   pgRect -60 -30 0 60
   pgRect -30 0 0 60
   pgRect 0 30 0 60
   pgRect 30 60 0 60
   pgSci 1
   }

#####################################################################
#Compute change in astigmatism with shift in focal plane.
#I would like to fetch surfaces of primary and focal plane using surface
#name, but alas not all designs have the requisite info yet.

proc astigShift {hndl primary focal filter} {
   set xmm [showFocal $hndl $filter xrad]
   if {$xmm < 0} {
	set ymm [showFocal $hndl $filter yrad]
   } else {
	set ymm 0.
	}
   set copy [opticNew]
   opticCopy $hndl $copy
   incrSurf $copy $focal z 2
   setColorFlag $copy $filter 1
   setSurfFlag $copy $primary z 1
   verbose 0
   global _optic
   set _optic(lscale) 0
   set _optic(fractmax) 1
   set _optic(niter) 1
   catch {unset _optic(flagcache)}
   clearSpot
   addSpot 0 0
   fastLstsq $copy
   set astig1 [zernikeShow $hndl $xmm $ymm $filter 2 2]
   set astig2 [zernikeShow $copy $xmm $ymm $filter 2 2]
   set a1 [expr sqrt(pow([lindex $astig1 0],2) + pow([lindex $astig1 1],2))]
   set a2 [expr sqrt(pow([lindex $astig2 0],2) + pow([lindex $astig2 1],2))]
   set slope [expr ($a2-$a1)/2.]
   opticDel $copy
   return $slope
   }

##################################################################
#For Jarvis and Jain, create a file with whiskers at finely tabulated
#intervals.

proc  jarvis {hndl ifil} {
   set xrad [showFocal $hndl $ifil xrad]
   set yrad [showFocal $hndl $ifil yrad]
   set scale [expr (.01*2.*$xrad/.1)]
   set nmax 100
   set name [showDesign $hndl]
   set file $name.txt
   set fid [open $file w]
   set format "%10.2f %10.2f %12.3f %12.3f %12.3f"
   set sformat "%10s %10s %12s %12s %12s"
   puts $fid [format $sformat xmm ymm e1(arcsec^2) e2(arcsec^2) \
	dilute(arcsec^2)]
   loop i 0 $nmax {
	set xmm [expr $xrad*($i-($nmax-1.)/2.)/($nmax/2.)]
	loop j 0 $nmax {
	   set ymm [expr $yrad*($j-($nmax-1.)/2.)/($nmax/2.)]
	   if {$xrad > 0 && $yrad > 0} {
		if {pow($xmm/$xrad,2) + pow($ymm/$yrad,2) > 1.} continue
		}
	   set list [polarize $hndl $xmm $ymm $ifil]
	   set e1 [lindex $list 0]
	   set e2 [lindex $list 1]
	   set dilute [lindex $list 2]
	   puts $fid [format $format $xmm $ymm $e1 $e2 $dilute]
	   }
	}
   close $fid
   return $file   }

#######################################################################
#Simulate a tilt of the prime focus corrector.  I do this by tilting
#and translating the primary mirror so that the pivot point is someplace
#close to the focal plane.
#
#Input is OPTIC, ifil, dz (positive) from focal plane towards primary
#I have hardwired which surfaces are the primary (2.0x) and focal plane (15)
#
#tilt is in arcsec
#The tilt and translation are all absolute - if I want to translate again,
#I do it after calling primTilt

proc primTilt {hndl ifil dz tilt} {
   set psurf 2.0[expr ($ifil+1)/2]

   set z2 [showSurf $hndl 15 z]
   set z1 [showSurf $hndl $psurf z]

#pivot is positive measured up from vertex of primary
   set pivot [expr $z1-$z2-$dz]
   set tilt [expr $tilt/(3600.*57.3)]

#If tilt is positive, I apply positive tilt and positive translation to primary
   set xmm [expr $tilt*$pivot]
echo Translation [format %.3f $xmm] mm
   setSurf $hndl $psurf x $xmm
   setSurf $hndl $psurf theta $tilt
   return
   }

