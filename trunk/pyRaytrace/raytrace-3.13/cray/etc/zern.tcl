#Set a Zernike term for a surface.  Input is n (order of radial polynomial)
#and m (order of angular polynomial).  Must have n even, n >= m.
#n an m must have same parity.
#For now, I will just set the cosine part
#targ will be a target rms spot size.  I will "rescale" to give amplitude
#of coefficient.
#
#I now provide an argument to select cosine or sine.  Default is cosine
#
# zset h0 surfid n m val [0|1]	Set a single Zernike coefficient on a surface
#				0 = cosine, 1 = sine
# zshow h0 surfid n m [0|1]	Show Zernike coeff
#				0 = cosine, 1 = sine
# zslist h0 surfid		Call zernikeList for the surface

# nzset h0 surfid n m val [0|1]	zset for refractive gradient zernike
# nzshow h0 surfid n m [0|1]	zshow for refractive gradient zernike
# nzlist h0 surfid		Call zernikeList for the refractive gradient
# zeval h0 surfid norder val	Test of scaling relations
# zpower h0 surfid norder val filter	Add a power spectrum of zernike's
#				up to an including norder and compute impact
#				on image size.  Image is at edge of field.
#				Must have run a least squares first to
#				set up compensator (refocus to get rid of sph.
#				aberr.)  Some randomization is done.

# zloop h0 surfid norder val filter	Run zpower 10 times and average
#				results.
# zrms h0 surfid n m		Test to make sure rms height matches theory.
#				Amplitude defaults to 1.
# zreg reg zern			Fill a region with a map computed from zernikes
# zscale zern scale		Scale coefficients of a zernike
# zernikeWrite zern file [reg]	Write out an ascii file of zernike terms
#				reg is region used for fitting zerikes
# zernikeRead file		Read back the ascii file
# zsetfromz optic surfid zern	Set zernike terms on a surface.

#NOTE ON COEFFIENT CONVENTION
#When analyzing wavefronts, my coefficients are in units of wavelength (waves
#of error).  For mirror surfaces, my coefficients are in mm.  When I run
#zernikeList, the column headings are incorrect for mirror surfaces (where
#it says waves, in actuality it is mm).

proc zset {optic surfid n m val {cs 0}} {
   set isurf [surfIndex $optic $surfid]
   if {$n < $m} {
	error "Radial order $n must be >= angular order $m"	
	}
   if {(($n-$m)/2)*2 != $n-$m} {
	error "n $n and m $m must have same parity"
	}
   if {$m == 0 && $cs > 0} return
   set N [expr $n+$m]
   set k [expr ($n-$m)/2]
   set j [expr ($N/2)*($N/2) + 2*$k + $cs]
   set zern [exprGet $optic.optic<$isurf>->zern]
   if {[regexp nil $zern]} {
	echo Initializing zernike structure
	opticZernikeAdd $optic $surfid $N
	}
   set ncoeff [exprGet $optic.optic<$isurf>->zern->ncoeff]
   if {$j >= $ncoeff} {
	echo Extending order of zernike structure
	opticZernikeAdd $optic $surfid $N
	}
   handleSet $optic.optic<$isurf>->zern->coeff<$j> $val
   return
   }

################################################################
proc zshow {optic surfid n m cs} {
   set isurf [surfIndex $optic $surfid]
   if {$n < $m} {
	error "Radial order $n must be >= angular order $m"	
	}
   if {(($n-$m)/2)*2 != $n-$m} {
	error "n $n and m $m must have same parity"
	}
   set N [expr $n+$m]
   set k [expr ($n-$m)/2]
   set j [expr ($N/2)*($N/2) + 2*$k + $cs]
   set zern [exprGet $optic.optic<$isurf>->zern]
   if {[regexp nil $zern]} {
	error "No zernike structure present"
	}
   set ncoeff [exprGet $optic.optic<$isurf>->zern->ncoeff]
   if {$j >= $ncoeff} {
	error "Specified term is out of bounds"
	}
   return [exprGet $optic.optic<$isurf>->zern->coeff<$j>]
   }

#######################################################################
#Same as above but for refractive index gradients.  These are are in units
#of ppm and are multipliers of the constant index.  (This helps avoid
#worrying about sign conventions).

proc nzset {optic surfid n m val {cs 0}} {
   set isurf [surfIndex $optic $surfid]
   if {$n < $m} {
	error "Radial order $n must be >= angular order $m"	
	}
   if {(($n-$m)/2)*2 != $n-$m} {
	error "n $n and m $m must have same parity"
	}
   if {$m == 0 && $cs > 0} return
   set N [expr $n+$m]
   set k [expr ($n-$m)/2]
   set j [expr ($N/2)*($N/2) + 2*$k + $cs]
   set zern [exprGet $optic.optic<$isurf>->nzern]
   if {[regexp nil $zern]} {
	echo Initializing zernike structure
	indexZernikeAdd $optic $surfid $N
	}
   set ncoeff [exprGet $optic.optic<$isurf>->nzern->ncoeff]
   if {$j >= $ncoeff} {
	echo Extending order of zernike structure
	indexZernikeAdd $optic $surfid $N
	}
   handleSet $optic.optic<$isurf>->nzern->coeff<$j> $val
   return
   }

################################################################
proc nzshow {optic surfid n m cs} {
   set isurf [surfIndex $optic $surfid]
   if {$n < $m} {
	error "Radial order $n must be >= angular order $m"	
	}
   if {(($n-$m)/2)*2 != $n-$m} {
	error "n $n and m $m must have same parity"
	}
   set N [expr $n+$m]
   set k [expr ($n-$m)/2]
   set j [expr ($N/2)*($N/2) + 2*$k + $cs]
   set zern [exprGet $optic.optic<$isurf>->nzern]
   if {[regexp nil $zern]} {
	error "No zernike structure present"
	}
   set ncoeff [exprGet $optic.optic<$isurf>->nzern->ncoeff]
   if {$j >= $ncoeff} {
	error "Specified term is out of bounds"
	}
   return [exprGet $optic.optic<$isurf>->nzern->coeff<$j>]
   }

#################################################################
#Do a test of scaling
#nmax is even

proc zeval {optic surfid nmax val} {
   echo [format "%2s %2s %2s %6s %6s %6s" N n m FWHM target ratio]
   for {set i 4} {$i <= $nmax} {incr i 2} {
	set power($i) 0.
	set tpower($i) 0.
	for {set m 0} {$m <= [expr $i/2]} {incr m} {
	   set n [expr $i-$m]
	   opticZernikeDel $optic $surfid
	   opticZernikeAdd $optic $surfid $nmax
	   zset $optic $surfid $n $m $val

#Target
	   set target [expr .1*($n+$m-1.)*($n-$m+1.)/$n]
	   if {$m > 0} {set target [expr $target/1.414]}
	   rtrace $optic 0 0 1 0
	   set fwhm [exprGet $optic.diagram->fwhmx]
	   set ratio [expr $fwhm/$target]
	   echo [format "%2d %2d %2d %.4f %.4f %.2f" $i $n $m $fwhm $target \
		$ratio]
	   set power($i) [expr $power($i) + $fwhm*$fwhm]
	   set tpower($i) [expr $tpower($i) + $target*$target]

#Double weight for m > 0 to allow for sine terms.
	   if {$m > 0} {
		set power($i) [expr $power($i) + $fwhm*$fwhm]
		set tpower($i) [expr $tpower($i) + $target*$target]
		}
	   }
	set power($i) [expr sqrt($power($i))]
	set tpower($i) [expr sqrt($tpower($i))]
	}
   foreach i [lsort -integer [array names power]] {
	echo Order $i total FWHM [format %.2f $power($i)] Target \
	   [format %.2f $tpower($i)]
	}
   return
   }

	
######################################################################
#Do a test of power for a given norder.  First, set all terms one by one
#and sum.  Then, set all simultaneously.  See if things sum OK.
#nmax is even
#Return 1-D FWHM in mm.

proc zpower {optic surfid norder amplitude {filter 1}} {

   set xmm [showFocal $optic $filter xrad]
   set amplitude [expr 1.*$amplitude]

#Clear out any old coefficients.
#Set all new coefficients to same value and see what pops out.
   opticZernikeDel $optic $surfid

#Refocus - must have set flags outside.
   fastLstsq $optic

#Get RMS FWHM before adding zernike's
   rtrace $optic $xmm 0 $filter 0
   set fwhmx [exprGet $optic.diagram->fwhmx]
   set fwhmy [exprGet $optic.diagram->fwhmy]
   set fwhm0 [expr sqrt(($fwhmx*$fwhmx + $fwhmy*$fwhmy)/2.)]

   opticZernikeAdd $optic $surfid $norder

#Now add zernike distortions.
#RMS Z height is analytically calculable.
   set norm 0.
   for {set iorder 4} {$iorder <= $norder} {incr iorder 2} {

	for {set m 0} {$m <= [expr $iorder/2]} {incr m} {
	   set n [expr $iorder-$m]

#Norm already included contribution of sine term for m>0.
	   set norm [expr $norm + 1./($n+1.)]

#Amplitude at order i.  If I scale like 1/iorder, the rms amplitude is capped,
#but the rms slope errors increase like sqrt(norder).
#
	   set val [expr $amplitude/$iorder]

#Apply random phase if m > 0
	   if {$m == 0} {
		set sign [random 2]
		if {$sign == 0} {
		   zset $optic $surfid $n $m [expr $val]
		} else {
		   zset $optic $surfid $n $m [expr -1.*$val]
		   }
	   } else {
		set phase [expr 2.*3.1416*[random 100]/100.]
		zset $optic $surfid $n $m [expr $val*cos($phase)] 0
		zset $optic $surfid $n $m [expr $val*sin($phase)] 1
		}
	   }
	}

#Refocus - must have set flags outside.
   fastLstsq $optic

#Evaluate PSF of distorted design.
   rtrace $optic $xmm 0 $filter 0
#   spotPlot $optic $xmm 0 $filter
   set fwhm 0.
   set fwhmx [exprGet $optic.diagram->fwhmx]
   set fwhmy [exprGet $optic.diagram->fwhmy]
   set fwhm [expr sqrt(($fwhmx*$fwhmx + $fwhmy*$fwhmy)/2.)]
   if {$fwhm > $fwhm0} {
	set fwhm [expr sqrt($fwhm*$fwhm - $fwhm0*$fwhm0)]
   } else {
	set fwhm 0
	}

#Analyze distortion pattern in detail.
#RMS z height and slope errors
   set rms 0.
   set slope 0.
   set zavg 0.
   set errxavg 0.
   set erryavg 0.
   set nfract [exprGet $optic.diagram->nfract]
   set isurf [surfIndex $optic $surfid]
   set outstop [expr abs([showSurf $optic $surfid outstop])];
   set minslope 999.
   set maxslope -999.
   loop i 0 $nfract {
	set xf [exprGet $optic.diagram->xfract<$i>]
	set yf [exprGet $optic.diagram->yfract<$i>]
	set xmm [expr $xf*$outstop]
	set ymm [expr $yf*$outstop]
	set list [zernikeSurfErr *$optic.optic<$isurf>->zern $xmm $ymm \
	   $outstop]
	set err [lindex $list 0]
	set errx [lindex $list 1]
	set erry [lindex $list 2]
	set rms [expr $rms + $err*$err]
	set slope [expr $slope + $errx*$errx + $erry*$erry]
	set minslope [min $minslope $errx $erry]
	set maxslope [max $maxslope $errx $erry]
	set zavg [expr $zavg + $err]
	set errxavg [expr $errxavg + $errx]
	set erryavg [expr $erryavg + $erry]
	}

#RMS is measured z rms; norm is analytic prediction.
   set rms [expr sqrt($rms/$nfract)]
   set zavg [expr ($zavg/$nfract)/.0006]
   set errxavg [expr ($errxavg/$nfract)*1.e6]
   set erryavg [expr ($erryavg/$nfract)*1.e6]
   set slope [expr sqrt($slope/$nfract)*1.e6]
   set pp [expr ($maxslope - $minslope)*1.e6]
   set norm [expr $val*sqrt($norm)]
#   echo rms z height is [format %.4f $rms], should be [format %.4f $norm]
   set waves [expr $rms/.0006]

#Convert FWHM to D80 in arcsec
   set scale [showFocal $optic $filter scale]
   set D80 [expr $fwhm*$scale*1.53]

#RMS surface error in microns
   set microns [expr $rms*1000.]
   echo Z rms (microns) [format %.3f $microns]
   echo Z rms (waves) [format %.2f $waves] D80 (arcsec) [format %.3f $D80] \
	Slope: RMS [format %.0f $slope] P-P [format %.0f \
	$pp] microrad

#Means are typpically quite small, so don't print.
#   echo Mean: Z [format %.2f $zavg] Slope-x [format %.0f $errxavg] Slope-y \
#	[format %.0f $erryavg]

#Cleanup.
#   mapZernike3d *$optic.optic<$isurf>->zern

#Metric: rms radius per wavefront error
   set rmsrad [expr $fwhm/2.35*1.414*1.e3]
   set xmm [showFocal $optic $filter xrad]
   set waveErr [waveFrontRms $optic $xmm 0 $filter]
   set metric [expr $rmsrad/$waveErr]
   echo RMS microns per wavefront error: [format %.1f $metric]

#Wave front rms.
   opticZernikeDel $optic $surfid
   return $D80
   }

######################################################################
#Loop 10 times through zpower and get average of D80
proc zloop {optic surfid norder amplitude {filter 1}} {
   set avg 0
   loop i 0 10 {
	set avg [expr $avg + [zpower $optic $surfid $norder $amplitude \
	   $filter]]
	}
   set avg [expr $avg/10.]
   echo Average D80 = [format %.3f $avg]
   return
   }

######################################################################
#Compute rms z position as a test.

proc zrms {optic surfid n m} {
   set NORDER [expr $m+$n]
   opticZernikeAdd $optic $surfid $NORDER
   zset $optic $surfid $n $m 1
   set rms 0.

#Use the current ray pattern
   set nfract [exprGet $optic.diagram->nfract]
   set isurf [surfIndex $optic $surfid]
   loop i 0 $nfract {
	set xf [exprGet $optic.diagram->xfract<$i>]
	set yf [exprGet $optic.diagram->yfract<$i>]
	set err [zernikeWaveErr *$optic.optic<$isurf>->zern $xf $yf]
	set rms [expr $rms + $err*$err]
	}
   set rms [expr sqrt($rms/$nfract)]

#What should it be?
   if {$m == 0} {
	set norm [expr 1./sqrt($n+1.)]
   } else {
	set norm [expr 1./sqrt(2.*($n+1.))]
	}
   echo rms [format %.2f $rms] should be [format %.2f $norm]
   return
   }

####################################################################
#List zernikes for a given surface.
#NOTE - there is a "proc zlist" in zernike.tcl which is intended, I think,
#for test pruposes only.  It is not loaded by default.
#Argh! there is also a zlist in "zlist.tcl"
proc zslist {hndl surf} {
   set isurf [surfIndex $hndl $surf]
   set zernike [exprGet $hndl.optic<$isurf>->zern]
   if {[regexp nil $zernike]} return
   zernikeList (*${hndl}.optic<$isurf>->zern)
   return
   }

####################################################################
#List zernikes for a given refractive gradient.
#NOTE - there is a "proc zlist" in zernike.tcl which is intended, I think,
#for test pruposes only.  It is not loaded by default.
proc nzlist {hndl surf} {
   set isurf [surfIndex $hndl $surf]
   set zernike [exprGet $hndl.optic<$isurf>->nzern]
   if {[regexp nil $zernike]} return
   zernikeList (*${hndl}.optic<$isurf>->nzern)
   return
   }

########################################################################
#Delete a zernike for a surface.  This is just a wrapper.

proc zdel {hndl surfid} {
   opticZernikeDel $hndl $surfid
   return
   }

########################################################################
#Delete a zernike for a refractive gradient.  This is just a wrapper.

proc nzdel {hndl surfid} {
   indexZernikeDel $hndl $surfid
   return
   }

#########################################################################
#Delete a generic zernike structure - just a wrapper

proc zernikeDel {zern} {
   zernikeZero $zern
   genericDel $zern
   return
   }

#########################################################################
#Set coeffs of a zernike structure directly - otherwise like zset and zshow
proc zzset {zern n m val {cs 0}} {
   if {$n < $m} {
	error "Radial order $n must be >= angular order $m"	
	}
   if {(($n-$m)/2)*2 != $n-$m} {
	error "n $n and m $m must have same parity"
	}
   if {$m == 0 && $cs > 0} return
   set N [expr $n+$m]
   set k [expr ($n-$m)/2]
   set j [expr ($N/2)*($N/2) + 2*$k + $cs]
   set ncoeff [exprGet $zern.ncoeff]
   if {$j >= $ncoeff} {
	echo Extending order of zernike structure
	opticZernikeAdd $optic $surfid $N
	}
   handleSet $zern.coeff<$j> $val
   return
   }

################################################################
proc zzshow {zern n m cs} {
   if {$n < $m} {
	error "Radial order $n must be >= angular order $m"	
	}
   if {(($n-$m)/2)*2 != $n-$m} {
	error "n $n and m $m must have same parity"
	}
   set N [expr $n+$m]
   set k [expr ($n-$m)/2]
   set j [expr ($N/2)*($N/2) + 2*$k + $cs]
   set ncoeff [exprGet $zern.ncoeff]
   if {$j >= $ncoeff} {
	error "Specified term is out of bounds"
	}
   return [exprGet $zern.coeff<$j>]
   }

##########################################################################
#Fill a region with surface errors computed via a zernike structure.
#Region must exist.

proc zreg {reg zern} {
   set nrow [exprGet $reg.nrow]
   set ncol [exprGet $reg.ncol]

   loop r 0 $nrow {
	loop c 0 $ncol {
	   set yfract [expr ($r + 0.5 -$nrow/2.)/($nrow/2.)]
	   set xfract [expr ($c + 0.5 -$ncol/2.)/($ncol/2.)]
	   set rad2 [expr $xfract*$xfract + $yfract*$yfract]
	   if {$rad2 >= 1.} continue
	   set err [zernikeWaveErr $zern $xfract $yfract]
	   set n [expr $r*$ncol + $c]
	   handleSet $reg.pixels<$n> $err
	   }
	}
   return
   }

#########################################################################
#Scale coefficients of a zernike.  Note that other auxiliary info (such as
#tabulated data points for a fit, etc.) are NOT scaled.
proc zscale {zern scale} {
   set ncoeff [exprGet $zern.ncoeff]
   loop j 0 $ncoeff {
	set val [exprGet $zern.coeff<$j>]
	set val [expr $val*$scale]
	handleSet $zern.coeff<$j> $val
	}
   return
   }

#########################################################################
#Write out coefficients (and NOTHING ELSE) to a file
#I really need to write out some provenance info - like where did the
#zernikes come from.  Normally I am deriving them from a fit to a region.
#I can now write comments at the start of the file and not break compatibility.

proc zernikeWrite {zern file {reg ""}} {
   set norder [exprGet $zern.norder]
   set ncoeff [exprGet $zern.ncoeff]
   set fid [open $file w]

#If we supply a region - assume that the zernikes are a fit to a wavefront
#map
   if {$reg != ""} {
	set name [exprGet $reg.name]
	set nrow [exprGet $reg.nrow]
	set ncol [exprGet $reg.ncol]
	puts $fid "#name $name"
	puts $fid "#nrow $nrow"
	puts $fid "#ncol $ncol"
	}
   puts $fid "norder $norder"
   puts $fid "ncoeff $ncoeff"
   loop j 0 $ncoeff {
	puts $fid "$j [exprGet $zern.coeff<$j>]"
	}
   close $fid
   return
   }

#########################################################################
#Read back coefficients (and NOTHING ELSE) from a file

proc zernikeRead {file} {
   set fid [open $file]

#Trim out leading comments
   while {1} {
	set line [string trim [gets $fid]]
	if {[string index $line 0] == "#"} continue
	break
	}
   if {[lindex $line 0] == "norder"} {
	set norder [lindex $line 1]
   } else {
	echo Gak! Expected norder but found $line
	set norder 0
	}
   set line [gets $fid]
   if {[lindex $line 0] == "ncoeff"} {
	set ncoeff [lindex $line 1]
   } else {
	echo Gak! Expected ncoeff but found $line
	set ncoeff 0
	}
   set zern [genericNew ZERNIKE]
   zernikeInit $zern $norder

   loop j 0 $ncoeff {
	set line [gets $fid]
	if {[lindex $line 0] != $j} {
	   echo Bad line: $line
	   break
	   }
	handleSet $zern.coeff<$j> [lindex $line 1]
	}
   close $fid
   return $zern
   }

########################################################################
#Given a Zerike array, assign coefficients to a zernike structure attached
#to a given optical surface.

proc zsetfromz {optic surfid zern} {

   set isurf [surfIndex $optic $surfid]
   set norder [exprGet $zern.norder]
   set ncoeff [exprGet $zern.ncoeff]
   opticZernikeAdd $optic $surfid $norder
   loop j 0 $ncoeff {

#Zero out tilt, focus terms
	if {$j <= 3} {
	   set val 0.
	} else {
	   set val [exprGet $zern.coeff<$j>]
	   }
	handleSet $optic.optic<$isurf>->zern->coeff<$j> $val
	}
   return
   }
