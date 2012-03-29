#Make 3d plots of various things.
#Procedures:
#   waveFrontMap3d optic xmm ymm icolor

#####################################################################
#The following routines is borrowed from plplot.
#
# Routine for initializing color map 1 in HLS space.
# Basic grayscale variation from half-dark (which makes more interesting
# looking plot compared to dark) to light.

proc cmap1_init_11 {} {

# Independent variable of control points.
   matrix i f 2 = {0., 1.}

# Hue ranges from blue (240 deg) to red (0 or 360 deg)
   matrix h f 2 = {240., 0.}

# Lightness and saturation are constant (values taken from C example).
   matrix l f 2 = {0.6, 0.6}
   matrix s f 2 = {0.8, 0.8}

# Integer flag array is zero (no interpolation along far-side of colour
# figure
   matrix rev i 2 = {0, 0}

# Number of cmap1 colours is 256 in this case. 
   plscmap1n 256

# Interpolate between control points to set up default cmap1.
   plscmap1l 0 2 i h l s rev
   return
   }

####################################################################
#Make a waveFront map and plot it.
#We remove mean slope from wavefront as well.

proc waveFrontMap3d {optic xmm ymm color} {

   set wave [genericNew WAVEFRONT]
   set NSTEP 10

#2 dimensions
   set waveMin 99
   set waveMax -99

   chiefRay $optic $wave $xmm $ymm $color

   matrix x f [expr 2*$NSTEP+1]
   matrix y f [expr 2*$NSTEP+1]
   matrix z f [expr 2*$NSTEP+1] [expr 2*$NSTEP+1]

   set mean 0
   set wavex 0
   set wavey 0
   set x2 0
   set y2 0
   set n 0
   set rms 0

   loop i 0 [expr 2*$NSTEP+1] {
	x $i = [expr ($i-$NSTEP)/double($NSTEP)]
	y $i = [expr ($i-$NSTEP)/double($NSTEP)]
	loop j 0 [expr 2*$NSTEP+1] {
	   set xfract [expr ($i-$NSTEP)/double($NSTEP)]
	   set yfract [expr ($j-$NSTEP)/double($NSTEP)]

#Error in waves
	   if {$xfract*$xfract+$yfract*$yfract <= .999} {
		set waveErr [offRay $optic $wave $xfract $yfract 0]
		set cache($i,$j) $waveErr
		set mean [expr $mean + $waveErr]
		set wavex [expr $wavex + $waveErr*$xfract]
		set wavey [expr $wavey + $waveErr*$yfract]
		set x2 [expr $x2 + $xfract*$xfract]
		set y2 [expr $y2 + $yfract*$yfract]
		incr n
	   } else {
		set waveErr 0.
		}
	   z $i $j = $waveErr
	   }
      }

   if {$n == 0} {
	error "No waves made it through system!"
	}
   set mean [expr $mean/$n]
   set a [expr $wavex/$x2]
   set b [expr $wavey/$y2]

   loop i 0 [expr 2*$NSTEP+1] {
	loop j 0 [expr 2*$NSTEP+1] {
	   set xfract [expr ($i-$NSTEP)/double($NSTEP)]
	   set yfract [expr ($j-$NSTEP)/double($NSTEP)]

#Error in waves
	   if {$xfract*$xfract+$yfract*$yfract <= .999} {
		set waveErr [expr $cache($i,$j) \
		   - $mean - $a*$xfract - $b*$yfract]
		set cache($i,$j) $waveErr
		set rms [expr $rms + $waveErr*$waveErr]
		set waveMin [min $waveMin $waveErr]
		set waveMax [max $waveMax $waveErr]
		}
	   }
      }

   set rms [expr sqrt($rms/$n)]

#Performance metrics in nanometers
   set waveLen [showWave $optic $color]
   set rmsnm [expr $rms*$waveLen*1.e3]
   set pp [expr $waveMax - $waveMin]
   set ppnm [expr $pp*$waveLen*1.e3]

   echo Wave Front Error (waves): [format %.3f $rms]
   loop i 0 [expr 2*$NSTEP+1] {
	x $i = [expr ($i-$NSTEP)/double($NSTEP)]
	y $i = [expr ($i-$NSTEP)/double($NSTEP)]
	loop j 0 [expr 2*$NSTEP+1] {
	   set xfract [expr ($i-$NSTEP)/double($NSTEP)]
	   set yfract [expr ($j-$NSTEP)/double($NSTEP)]

#Error in waves
	   if {$xfract*$xfract+$yfract*$yfract <= .999} {
		set waveErr $cache($i,$j)
	   } else {
		set waveErr $waveMin
		}
	   z $i $j = $waveErr
	   }
      }

   genericDel $wave
   set zmin $waveMin
   set zmax $waveMax
   if {$waveMax - $waveMin < 1} {
	set zmax [expr $zmin+1]
	}
   plotInit a
   pladv 0
   plcol0 15
   plvpor 0.0 1.0 0.0 0.9
   plwind -1.0 1.0 -1.0 1.5

   cmap1_init_11
   plw3d 1.0 1.0 1.2 -1.1 1.1 -1.1 1.1 $zmin $zmax 40. 60
   plbox3 "bnstu" "x axis" 0.0 0 \
	"bnstu" "y axis" 0.0 0 \
	"bcdmnstuv" "z axis" 0.0 4
   plcol0 2

   set MAG_COLOR 0x04
#   plmesh x y z 7
   plot3d x y z 6 0
   pgSci 1
   pgLabel "" "" "Wavefront Analysis: xmm = $xmm, ymm = $ymm"
   pgText -.9 -.75 "Peak-Peak: [format %.2f $pp] waves"
   pgText 0 -.75 "RMS: [format %.2f $rms] waves"
   pgText -.9 -.88 "             [format %.0f $ppnm] nm"
   pgText 0 -.88 "      [format %.0f $rmsnm] nm"

#Estimated Strehl ratio (Marechal relation)
   set strehl [expr exp(-pow(2.*3.1416*$rms,2))]
   echo Estimated Strehl ratio [format %.2f $strehl]
   return
   }

######################################################################
#Make a 3d mesh plot of a diffraction pattern.

proc psfMap3d {optic xmm ymm icolor} {
   global NPIX NPUPIL PIXMM SCALE
   if {[info exists NPUPIL]} {unset NPUPIL}
   if {[info exists PIXMM]} {unset PIXMM}
   if {[info exists SCALE]} {unset SCALE}
   set NPIX 32
   psfMap $optic $xmm $ymm $icolor
   set file fft$icolor.fit
   set reg [regReadFromFits $file]
   set pixmax 0
   set nrow [exprGet $reg.nrow]
   set ncol [exprGet $reg.ncol]
   matrix x f $nrow
   matrix y f $ncol
   matrix z f $nrow $ncol
   loop i 0 $nrow {
	x $i = $i
	}
   loop i 0 $ncol {
	y $i = $i
	}
   loop i 0 $nrow {
	loop j 0 $ncol {
	   set pix [expr sqrt([regPixGet $reg $i $j])]
	   z $i $j = $pix
	   set pixmax [max $pixmax $pix]
	   }
	}
   regDel $reg
   plotInit b
   pladv 0
   plcol0 15
   plvpor 0.0 1.0 0.0 0.9
   plwind  -1.0 1.0 -1 1.5

   cmap1_init_11
   plw3d 1.0 1.0 1.2 0 $nrow 0 $ncol 0 $pixmax 40 60
   plbox3 "bnstu" "x axis" 0.0 0 \
	"bnstu" "y axis" 0.0 0 \
	"bcdmnstuv" "z axis" 0.0 4
   plcol0 2

   set MAG_COLOR 0x04
#   plmesh x y z 7
   plot3d x y z 6 0
   pgSci 1
   pgLabel "" "" "PSF Analysis: xmm = $xmm, ymm = $ymm, wave = \
	[showFocal $optic $icolor wave]"
   return
   }

#################################################
#Fit a zernike to a wavefront.

proc zernikeWaveFit {optic xmm ymm color} {

   set zern [genericNew ZERNIKE]
   set thresh .005

#Find lowest order that meets threshold test.
   set NORDER 2
   while {1} {
        zernikeFit $optic $zern $NORDER $xmm $ymm $color
        set err [exprGet $zern.fitErr]
        if {$err < $thresh} break

#Sometimes we do want to fit to still higher order - e.g., evaluating errors
#on C4.  I will not change for now.
        if {$NORDER >= 14} {
	   echo NORDER > 14
	   break
	   }
        incr NORDER 2
        }
   return $zern
   }

####################################################################
#Make a zernike fit and return values for the specified order.
#Return in microns.

proc zernikeShow {optic xmm ymm color n m} {
   if {$n < $m} {
	error "Radial order $n must be >= angular order $m"
	}
   if {(($n-$m)/2)*2 != $n-$m} {
	error "n $n and m $m must have same parity"
	}
   set N [expr $n+$m]
   set k [expr ($n-$m)/2]
   if {$m > 0} {
	set nterm 2
   } else {
	set nterm 1
	}
   set list ""
   set zern [zernikeWaveFit $optic $xmm $ymm $color]
   set wave [exprGet $zern.wavelen]
   loop cs 0 $nterm {
	set j [expr ($N/2)*($N/2) + 2*$k + $cs]
	set ncoeff [exprGet $zern.ncoeff]
	if {$j >= $ncoeff} {
	   lappend list 0
	} else {
	   set waves [exprGet $zern.coeff<$j>]
	   set microns [expr $waves*$wave]
	   lappend list $microns
	   }
	}
   zernikeZero $zern
   genericDel $zern
   return $list
   }

####################################################################
#Make a waveFront map and plot it.  Use Zernike's, and wipe out the
#first 3 coefficients to get rid of constant and linear slopes (which are
#inconsequential for system performance)

proc zernikeMap3d {optic xmm ymm color} {

   set NSTEP 10

#Setting up zernike array and making fit is now broken out in a separate
#routine
   set zern [zernikeWaveFit $optic $xmm $ymm $color]

#Wipe out the constant and cos(theta), sin(theta) terms.  These are 0, 1, 2 in
#vec array
   handleSet $zern.coeff<0> 0
   handleSet $zern.coeff<1> 0
   handleSet $zern.coeff<2> 0

#2 dimensions
   set waveMin 99
   set waveMax -99

   matrix x f [expr 2*$NSTEP+1]
   matrix y f [expr 2*$NSTEP+1]
   matrix z f [expr 2*$NSTEP+1] [expr 2*$NSTEP+1]
   set mean 0
   set rms 0
   set n 0
   loop i 0 [expr 2*$NSTEP+1] {
	x $i = [expr ($i-$NSTEP)/double($NSTEP)]
	y $i = [expr ($i-$NSTEP)/double($NSTEP)]
	loop j 0 [expr 2*$NSTEP+1] {
	   set xfract [expr ($i-$NSTEP)/double($NSTEP)]
	   set yfract [expr ($j-$NSTEP)/double($NSTEP)]
#Error in waves
	   if {$xfract*$xfract+$yfract*$yfract <= .999} {
		set waveErr [zernikeWaveErr $zern $xfract $yfract]
		set mean [expr $mean + $waveErr]
		set rms [expr $rms + $waveErr*$waveErr]
		incr n
		set waveMin [min $waveMin $waveErr]
		set waveMax [max $waveMax $waveErr]
		}
	   }
      }
   if {$n == 0} {
	error "No waves made it through system!"
	}
   set mean [expr $mean/$n]

   set mean 0
   set rms [expr sqrt($rms/$n - $mean*$mean)]
   set waveMin [expr $waveMin - $mean]
   set waveMax [expr $waveMax - $mean]

#Performance metrics in nanometers
   set waveLen [showWave $optic $color]
   set rmsnm [expr $rms*$waveLen*1.e3]
   set pp [expr $waveMax - $waveMin]
   set ppnm [expr $pp*$waveLen*1.e3]

#Print out waveFront error from zernike calculation.
#Raw errror is computed using actual wavefront; fitted error is from
#zernike approximation to wavefront.
#   echo Raw Wave Front Error (waves): [format %.3f [exprGet $zern.rmsErr]]
   echo Fit Wave Front Error (waves): [format %.3f $rms]

   loop i 0 [expr 2*$NSTEP+1] {
	x $i = [expr ($i-$NSTEP)/double($NSTEP)]
	y $i = [expr ($i-$NSTEP)/double($NSTEP)]
	loop j 0 [expr 2*$NSTEP+1] {
	   set xfract [expr ($i-$NSTEP)/double($NSTEP)]
	   set yfract [expr ($j-$NSTEP)/double($NSTEP)]

#Error in waves
	   if {$xfract*$xfract+$yfract*$yfract <= .999} {
		set waveErr [expr [zernikeWaveErr $zern $xfract \
		$yfract] - $mean]
	   } else {
		set waveErr $waveMin
		}
	   z $i $j = $waveErr
	   }
      }
   zernikeList $zern
   zernikeZero $zern
   genericDel $zern
   set zmin $waveMin
   set zmax $waveMax
   if {$waveMax - $waveMin < 1} {
	set zmax [expr $zmin+1]
	}
   plotInit a
   pladv 0
   plcol0 15
   plvpor 0.0 1.0 0.0 0.9
   plwind -1.0 1.0 -1.0 1.5

   cmap1_init_11
   plw3d 1.0 1.0 1.2 -1.1 1.1 -1.1 1.1 $zmin $zmax 40. 60
   plbox3 "bnstu" "x axis" 0.0 0 \
	"bnstu" "y axis" 0.0 0 \
	"bcdmnstuv" "z axis" 0.0 4
   plcol0 2

   set MAG_COLOR 0x04
#   plmesh x y z 7
   plot3d x y z 6 0
   pgSci 1
   pgLabel "" "" "Wavefront Analysis: xmm = $xmm, ymm = $ymm"
   pgText -.9 -.75 "Peak-Peak: [format %.2f $pp] waves"
   pgText 0 -.75 "RMS: [format %.2f $rms] waves"
   pgText -.9 -.88 "             [format %.0f $ppnm] nm"
   pgText 0 -.88 "      [format %.0f $rmsnm] nm"

#Estimated Strehl ratio (Marechal relation)
   set strehl [expr exp(-pow(2.*3.1416*$rms,2))]
   echo Estimated Strehl ratio [format %.2f $strehl]
   return
   }

####################################################################
#Compute rms waveFront error and use as a figure of merit.
#Use the zernike routines, since they compute W.F.E. using a dense grid of
#rays, normally denser than what I input through rayPattern.
#Note that this routine does NOT obey stops.

proc waveFrontRms {optic xmm ymm color} {

   set zern [genericNew ZERNIKE]

#Find lowest order that meets threshold test.
   set NORDER 2
   zernikeFit $optic $zern $NORDER $xmm $ymm $color

   zernikeZero $zern
   set rmsErr [exprGet $zern.rmsErr]
   genericDel $zern
   return $rmsErr
   }

####################################################################
#Map an existing Zernike.

proc mapZernike3d {zern} {
   set NSTEP 10

#2 dimensions
   set waveMin 99
   set waveMax -99

   matrix x f [expr 2*$NSTEP+1]
   matrix y f [expr 2*$NSTEP+1]
   matrix z f [expr 2*$NSTEP+1] [expr 2*$NSTEP+1]
   set mean 0
   set rms 0
   set n 0
   loop i 0 [expr 2*$NSTEP+1] {
	x $i = [expr ($i-$NSTEP)/double($NSTEP)]
	y $i = [expr ($i-$NSTEP)/double($NSTEP)]
	loop j 0 [expr 2*$NSTEP+1] {
	   set xfract [expr ($i-$NSTEP)/double($NSTEP)]
	   set yfract [expr ($j-$NSTEP)/double($NSTEP)]
#Error in waves
	   if {$xfract*$xfract+$yfract*$yfract <= .999} {
		set waveErr [zernikeWaveErr $zern $xfract $yfract]
		set mean [expr $mean + $waveErr]
		set rms [expr $rms + $waveErr*$waveErr]
		incr n
		set waveMin [min $waveMin $waveErr]
		set waveMax [max $waveMax $waveErr]
		}
	   }
      }
   if {$n == 0} {
	error "No waves made it through system!"
	}
   set mean [expr $mean/$n]

   set mean 0
   set rms [expr sqrt($rms/$n - $mean*$mean)]
   set waveMin [expr $waveMin - $mean]
   set waveMax [expr $waveMax - $mean]

#Performance metrics in nanometers
   set waveLen .5
   set rmsnm [expr $rms*$waveLen*1.e3]
   set pp [expr $waveMax - $waveMin]
   set ppnm [expr $pp*$waveLen*1.e3]

#Print out waveFront error from zernike calculation.
#Raw errror is computed using actual wavefront; fitted error is from
#zernike approximation to wavefront.
#   echo Raw Wave Front Error (waves): [format %.3f [exprGet $zern.rmsErr]]
   echo Fit Wave Front Error (waves): [format %.3f $rms]

   loop i 0 [expr 2*$NSTEP+1] {
	x $i = [expr ($i-$NSTEP)/double($NSTEP)]
	y $i = [expr ($i-$NSTEP)/double($NSTEP)]
	loop j 0 [expr 2*$NSTEP+1] {
	   set xfract [expr ($i-$NSTEP)/double($NSTEP)]
	   set yfract [expr ($j-$NSTEP)/double($NSTEP)]

#Error in waves
	   if {$xfract*$xfract+$yfract*$yfract <= .999} {
		set waveErr [expr [zernikeWaveErr $zern $xfract \
		$yfract] - $mean]
	   } else {
		set waveErr $waveMin
		}
	   z $i $j = $waveErr
	   }
      }
   set zmin $waveMin
   set zmax $waveMax
   if {$waveMax - $waveMin < 1} {
	set zmax [expr $zmin+1]
	}
   plotInit a
   pladv 0
   plcol0 15
   plvpor 0.0 1.0 0.0 0.9
   plwind -1.0 1.0 -1.0 1.5

   cmap1_init_11
   plw3d 1.0 1.0 1.2 -1.1 1.1 -1.1 1.1 $zmin $zmax 40. 60
   plbox3 "bnstu" "x axis" 0.0 0 \
	"bnstu" "y axis" 0.0 0 \
	"bcdmnstuv" "z axis" 0.0 4
   plcol0 2

   set MAG_COLOR 0x04
#   plmesh x y z 7
   plot3d x y z 6 0
   pgSci 1
   pgText -.9 -.75 "Peak-Peak: [format %.2f $pp] waves"
   pgText 0 -.75 "RMS: [format %.2f $rms] waves"
   pgText -.9 -.88 "             [format %.0f $ppnm] nm"
   pgText 0 -.88 "      [format %.0f $rmsnm] nm"

#Estimated Strehl ratio (Marechal relation)
   set strehl [expr exp(-pow(2.*3.1416*$rms,2))]
   echo Estimated Strehl ratio [format %.2f $strehl]
   return
   }

####################################################################
#List the contents of a zernike structure along with some 
#helper information.

proc zernikeList {zern} {
   set NORDER [exprGet $zern.norder]
   set ncoeff [exprGet $zern.ncoeff]

   set M [expr $NORDER/2]
   set jlist ""
   for {set m 0} {$m <= $M} {incr m} {
	for {set n $m} {$n <= $NORDER-$m} {set n [expr $n+2]} {
	   set k [expr ($n-$m)/2]
	   set N [expr ($n+$m)]
	   set j [expr ($N/2)*($N/2) + 2*$k]
	   lappend jlist $j
	   set table($j,n) $n
	   set table($j,m) $m
	   if {$m > 0} {
		incr j
		lappend jlist $j
		set table($j,n) $n
		set table($j,m) $m
		}
	   }
	}

#Popular names
   set name(0) Piston
   set name(1) "Tilt about x axis"
   set name(2) "Tilt about y axis"
   set name(3) Defocus
   set name(4) "Astigmatism with axis at 0 or 90 deg"
   set name(5) "Astigmatism with axis at +/- 45 deg"
   set name(6) "Coma along x axis"
   set name(7) "Coma along y axis"
   set name(8) "Spherical aberration"
   set name(9) "Triangular astig., base on x axis"
   set name(10) "Triangular stig., base on y axis"
   set name(16) "Ashtray astig., crests on axis"
   set name(17) "Ashtray astig., nodes on axes"

   puts stdout [format "%4s %5s  %5s  %10s  %10s  %s" \
	Term Rad. Ang. Waves Microns Name]
   set format "%3d  %5d  %5d  %10.3f  %10.3f  %-s"
   set wave [exprGet $zern.wavelen]
   foreach j [lsort -integer $jlist] {
	if {[info exists name($j)]} {
	   set info $name($j)
	} else {
	   set info ""
	   }
	set waves [exprGet $zern.coeff<$j>]
	set microns [expr $waves*$wave]
	puts stdout [format $format $j $table($j,n) $table($j,m) \
	   $waves $microns $info]
	}
   return
   }

