#Using zernikeMap:
#	set PIXMM ...
#	set NPIX ...
#	set NPUPIL ...
#Top level procedures for calling FFT diffraction analysis routines.
#Note: These procedure names are not really informative, but I already used
#fftMap, etc in calls to C code, so I'm at a loss for anything more creative.

#First, a procedure to find the appropriate NORDER for a given spot
#position.  Input a threshold for rms wavefront error.
#Threshold of .002 seems pretty good

####################################################################
#Call psfMap, but just print out F, G for checking
proc psfCheck {optic xmm ymm colors} {
   psfMap $optic $xmm $ymm $colors 0
   return
   }
####################################################################
#Compute diffraction map using FFT method.

proc psfMap {optic xmm ymm colors {mode 1}} {
   global THRESH NPIX SCALE PIXMM NPUPIL SIGMA
   if {[info exists THRESH]} {
	set thresh $THRESH
   } else {
	set thresh .002
	}
   if {[info exists NPIX]} {
	set npix $NPIX
   } else {
	set npix 128
	}
   set G 2.

#Telescope params.
   opticInfo $optic [lindex $colors 0]

#Focal length
   set fl [expr abs([showFocal $optic [lindex $colors 0] fl])]

#Mirror diameter - this is twice the radius of the aperture stop.
   set diam [telDiam $optic]

#Final focal ratio
   set fr [expr $fl/$diam]
   echo Focal ratio: $fr

#What is the scale per pixel?
#Scale is set by the shortest wavelength
   set minwave 9999.
   foreach icolor $colors {
	set minwave [min $minwave [showFocal $optic $icolor wave]]
	}

#Arcsec per pixel.
   set scale [expr ($minwave*1.e-3/$diam)*(1./$G)*206265.]
   if {[info exists SCALE]} {
	set scale $SCALE
	}
   set pixmm [expr $scale*$fl/206265.]
   if {[info exists PIXMM]} {
	set pixmm $PIXMM
	set scale [expr $pixmm*206265./$fl]
	}
echo pixel size(mm) [format %.4f $pixmm]
   set field [expr $scale*$npix]
   set fieldmm [expr $pixmm*$npix]
   echo Field size = [format %.2f $field] arcsec = [format %.3f $fieldmm] mm.
   echo NPIX = $npix
   if {[info exists SIGMA]} {
	echo Convolving with Gaussian, sigma = $SIGMA
	}
   set maxxcen -1.e10
   set minxcen 1.e10
   set maxycen -1.e10
   set minycen 1.e10
   foreach icolor $colors {
	set wave [showFocal $optic $icolor wave]
	if {[info exists SIGMA]} {
	   set sig [expr $SIGMA/$pixmm]
	} else {
	   set sig 0
	   }
	set G [expr $wave*1.e-3*206265./($diam*$scale)]

#What is F?  Let me try inputting the number of points across the exit pupil.
#The number I get for F=1 is NPIX/G.
	set F [expr int($G+.5)]

	if {[info exists NPUPIL]} {
	   set F [expr int($NPUPIL*$G/$npix + 0.5)]
	   if {$F < 1} {set F 1}
	   }
	set npupil($icolor) [expr int($F*$npix/$G)]
	echo G = [format %.1f $G] F = $F NPUPIL = $npupil($icolor)
	if {$mode == 0} continue
	set reg($icolor) [fftMap $optic $xmm $ymm $npix $F $G $icolor]
	if {$sig != 0} {fftConvolve $reg($icolor) $sig}
	hdrInsWithDbl $reg($icolor).hdr STREHLNM \
	   [exprGet $reg($icolor).param] "Central Intensity of Perfect Image"

#Compute x,y of the chief ray
	ray $optic $xmm $ymm 0 0 $icolor 0
	set index [expr [exprGet $optic.diagram->np] - 1]
	set xref($icolor) [expr [exprGet $optic.diagram->xray<$index>]/$pixmm]
	set yref($icolor) [expr [exprGet $optic.diagram->yray<$index>]/$pixmm]
	set maxxcen [max $maxxcen $xref($icolor)]
	set minxcen [min $minxcen $xref($icolor)]
	set maxycen [max $maxycen $yref($icolor)]
	set minycen [min $minycen $yref($icolor)]
	}
   if {$mode == 0} return

#Subregion field size
   set icolor1 [lindex $colors 0]
   set ncol [expr int([exprGet $reg($icolor1).ncol] - ($maxxcen-$minxcen))]
   set nrow [expr int([exprGet $reg($icolor1).nrow] - ($maxycen-$minycen))]
   foreach icolor $colors {
	set xoff [expr $minxcen - $xref($icolor)]
	if {$xoff < 0} {set xoff 0}
	if {$xoff+$ncol > [exprGet $reg($icolor).ncol]} {
	   set xoff [expr [exprGet $reg($icolor).ncol]-$ncol]
	   }
	set yoff [expr $minycen - $xref($icolor)]
	if {$yoff < 0} {set yoff 0}
	if {$yoff+$nrow > [exprGet $reg($icolor).nrow]} {
	   set yoff [expr [exprGet $reg($icolor).nrow]-$nrow]
	   }
	set subreg($icolor) [subRegNew $reg($icolor) $nrow $ncol $yoff $xoff]
	set norm [hdrGetAsDbl $reg($icolor).hdr STREHLNM]
	hdrInsertLine $subreg($icolor).hdr 1 \
"COMMENT FFT-enhanced diffraction image from CRAY fftMap"
	global env
	hdrInsWithAscii $subreg($icolor).hdr VERSION \
	   [file tail $env(CRAY_DIR)] "CRAY version"
	hdrInsWithAscii $subreg($icolor).hdr DESIGN \
	   [lindex [exprGet $optic.name] 0]
	hdrInsWithInt $subreg($icolor).hdr FILTER $icolor "Filter index"
	hdrInsWithDbl $subreg($icolor).hdr WAVE [showFocal $optic \
	   $icolor wave] "Wavelength (microns)"
	hdrInsWithDbl $subreg($icolor).hdr SCALE $scale "Arcsec/pixel"
	hdrInsWithDbl $subreg($icolor).hdr PIXMM $pixmm "Pixel size (mm)"
	hdrInsWithInt $subreg($icolor).hdr NPIX $npix \
	   "Number of pixels in original image"

#npupil includes subsampling and is total number of samples across a diameter.
	hdrInsWithInt $subreg($icolor).hdr NPUPIL $npupil($icolor) \
	   "Number of samples across exit pupil"
	hdrInsWithDbl $subreg($icolor).hdr XMM $xmm \
	   "Target x position (mm) in focal plane"
	hdrInsWithDbl $subreg($icolor).hdr YMM $ymm \
	   "Target y position (mm) in focal plane"

#Strehl stuff
	hdrInsWithDbl $subreg($icolor).hdr STREHLNM $norm "Perfect image \
central intensity"
###	set pix [regPixGet $subreg($icolor) [expr $nrow/2] [expr $ncol/2]]
	set list [regStatsFind $subreg($icolor)]
	set pix [keylget list high]
	set strehlratio [expr $pix/$norm]
	hdrInsWithDbl $subreg($icolor).hdr STREHL $strehlratio \
	   "Strehl ratio"

#Total flux in region.  At the moment, this is hardwired in fftmap.c
	hdrInsWithDbl $subreg($icolor).hdr FLUX 1000. \
	   "Total design flux in region"
	echo Filter $icolor Strehl ratio [format %.2f $strehlratio]
	regWriteAsFits $subreg($icolor) fft$icolor.fit
	regDel $subreg($icolor)
	regDel $reg($icolor)
	echo fft$icolor.fit
	}
   return
   }

####################################################################
#Create PPM from fits files.

proc fileToPPM {file1 file2 file3 outfile} {
   set reg1 [regReadAsFits [regNew] $file1]
   set reg2 [regReadAsFits [regNew] $file2]
   set reg3 [regReadAsFits [regNew] $file3]
   regToPPM $reg1 $reg2 $reg3 $outfile
   regDel $reg1
   regDel $reg2
   regDel $reg3
   return
   }

####################################################################
#Compute diffraction map using direct computation, but with no wavefront
#errors.  Useful for computing large scale diffraction pattern.

proc directMap {optic xmm ymm colors zstop nstrut strutw {strutang 0}} {
   global THRESH NPIX SCALE PIXMM NPUPIL SIGMA
   if {[info exists NPIX]} {
	set npix $NPIX
   } else {
	set npix 128
	}

#Telescope params.
   opticInfo $optic [lindex $colors 0]

#Focal length
   set fl [expr abs([showFocal $optic [lindex $colors 0] fl])]

#Mirror diameter - this is twice the radius of the aperture stop.
   set diam [telDiam $optic]

#Radius
   set radius [expr $diam/2.]

#Final focal ratio
   set fr [expr $fl/$diam]
   echo Focal ratio: $fr

#Inner radius
   set finner [exprGet $optic.tel->finner]
   echo Fractional inner radius $finner

#Pixel size (default is 10 microns)
   if {[info exists PIXMM]} {
	set pixmm $PIXMM
   } else {
	set pixmm .01
	}
   set scale [showFocal $optic [lindex $colors 0] scale]
   echo pixel size(mm) [format %.4f $pixmm]
   set field [expr $scale*$npix*$pixmm]
   set fieldmm [expr $pixmm*$npix]
   echo Field size = [format %.2f $field] arcsec = [format %.3f $fieldmm] mm.
   echo NPIX = $npix
   if {[info exists SIGMA]} {
	echo Convolving with Gaussian, sigma = $SIGMA
	}
   foreach icolor $colors {
	if {[info exists SIGMA]} {
	   set sig [expr $SIGMA/$pixmm]
	} else {
	   set sig 0
	   }
	set reg($icolor) [diffractMap $optic $xmm $ymm $icolor $npix $pixmm \
	   $fr $radius $zstop $finner $nstrut $strutw $strutang]
	if {$sig != 0} {fftConvolve $reg($icolor) $sig}

#Compute x,y of the chief ray
	ray $optic $xmm $ymm 0 0 $icolor 0
	set index [expr [exprGet $optic.diagram->np] - 1]
	set xref($icolor) [expr [exprGet $optic.diagram->xray<$index>]/$pixmm]
	set yref($icolor) [expr [exprGet $optic.diagram->yray<$index>]/$pixmm]
	}

   foreach icolor $colors {

	global env
	hdrInsWithAscii $reg($icolor).hdr VERSION \
	   [file tail $env(CRAY_DIR)] "CRAY version"
	hdrInsWithAscii $reg($icolor).hdr DESIGN \
	   [lindex [exprGet $optic.name] 0]
	hdrInsWithInt $reg($icolor).hdr FILTER $icolor "Filter index"
	hdrInsWithDbl $reg($icolor).hdr WAVE [showFocal $optic \
	   $icolor wave] "Wavelength (microns)"
	hdrInsWithDbl $reg($icolor).hdr PIXMM $pixmm "Pixel size (mm)"

	hdrInsWithDbl $reg($icolor).hdr XMM $xmm \
	   "Target x position (mm) in focal plane"
	hdrInsWithDbl $reg($icolor).hdr YMM $ymm \
	   "Target y position (mm) in focal plane"

#Total flux in region.  At the moment, this is hardwired in diffract.c
	hdrInsWithDbl $reg($icolor).hdr FLUX 1000. \
	   "Total design flux in region"
	regWriteAsFits $reg($icolor) fft$icolor.fit
	regDel $reg($icolor)
	echo fft$icolor.fit
	}
   return
   }

########################################################
#Modulation transfer function.
#Assume that I have already created a region with the PSF.  Pass in to
#here.  Output is a 1-D vector with the MTF

proc psfMtf {file} {
   set reg [regReadFromFits $file]
   set nrow [exprGet $reg.nrow]
   set ncol [exprGet $reg.ncol]

   fftMtf $reg

   set pixscale [hdrGetAsDbl $reg.hdr SCALE]

#Units are cycles per (pixscale*npix) arcsec
#(or something like that)

#Create an azimuthally-averaged profile.
#This is a real hack.
   set nstep [expr $nrow/2]
   loop i 0 $nstep {
	set npix($i) 0
	set pix($i) 0.
	}
   loop i 0 [expr $nrow/2] {
	loop j 0 [expr $ncol/2] {
	   set p [regPixGet $reg $i $j]
	   set irad [expr round(sqrt($i*$i*1. + $j*$j*1.))]
	   if {$irad >= $nstep} continue
	   incr npix($irad)
	   set pix($irad) [expr $pix($irad) + $p]
	   }
	}
   loop i 0 $nstep {
	if {$npix($i) > 0} {
	   set pix($i) [expr $pix($i)*1./$npix($i)]
	   }
	}
   set scale $pix(0)
   set xlist ""
   set ylist ""

#This gives me cycles per arcsec.
   set kstep [expr 1./($pixscale*$nrow)]

#Schroeder uses a normalized spatial frequency, where the units are
#D/lambda arcsec^(-1)
#This gives me cycles per pixel.
#For good sampling, we want cycles per pixel to go to 0 at 0.5 or lower.
#   set kstep [expr 1./($nrow)]

   loop i 0 $nstep {
	set pix($i) [expr $pix($i)/$scale]
	set x [expr $i*$kstep]
	lappend xlist $x
	lappend ylist $pix($i)
	}
   set xmax [expr $nstep*$kstep]
   pgEnv 0 $xmax 0 1.05 0 0
   pgLine $xlist $ylist
   pgPoint $xlist $ylist 3
   pgLabel "Cycles per arcsec" "MTF" $file
   regDel $reg
   return
   }
