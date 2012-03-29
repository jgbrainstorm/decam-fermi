#Compute and list the Seidel aberrations

proc opticSeidel {optic xmm ymm icolor {icolor2 0}} {
   seidel $optic $xmm $ymm $icolor $icolor2
   set format "%6s %9s %9s %9s %9s"
   set format2 " %9s %9s"
   puts -nonewline stdout [format $format SurfId TAS TSC TSA PETZ]
   if {$icolor2 != 0} {puts -nonewline stdout [format $format2 TAC TLC]}
   puts stdout " (microns)"
   set format "%6s %9.0f %9.0f %9.0f %9.0f"
   set format2 " %9.0f %9.0f"
   set tassum 0.
   set tscsum 0.
   set tsasum 0.
   set tacsum 0.
   set tlcsum 0.
   set hsum 0.
   set r [expr sqrt($xmm*$xmm + $ymm*$ymm)]
   loop i 1 [expr [exprGet $optic.diagram->np] - 1] {
	set isurf [exprGet $optic.diagram->indx<$i>]
	set surfid [surfId $optic $isurf]
	set tas [expr [exprGet $optic.optic<$isurf>->tas]*1.e3]
	set tsc [expr [exprGet $optic.optic<$isurf>->tsc]*1.e3]
	set tsa [expr [exprGet $optic.optic<$isurf>->tsa]*1.e3]
	set petz [exprGet $optic.optic<$isurf>->petz]
	set tac [expr [exprGet $optic.optic<$isurf>->tac]*1.e3]
	set tlc [expr [exprGet $optic.optic<$isurf>->tlc]*1.e3]

#Express field curvature as height of focal plane at edge.
	set height [expr .5*$petz*$r*$r*1.e3]
	set tassum [expr $tassum + $tas]
	set tscsum [expr $tscsum + $tsc]
	set tsasum [expr $tsasum + $tsa]
	set tacsum [expr $tacsum + $tac]
	set tlcsum [expr $tlcsum + $tlc]
	set hsum [expr $hsum + $height]


#Typical units are microns.
	puts -nonewline stdout [format $format $surfid $tas $tsc $tsa $height]
	if {$icolor2 != 0} {
	   puts -nonewline stdout [format $format2 $tac $tlc]
	   }
	puts stdout ""
	}
   puts -nonewline stdout [format $format Sum $tassum $tscsum $tsasum $hsum]
   if {$icolor2 != 0} {
	puts -nonewline stdout [format $format2 $tacsum $tlcsum]
	}
   puts stdout ""
   return
   }

############################################################################
#Same as opticSeidel, but groups surfaces together in lenses.
#For now, I will take a short-cut and group all surfaces with the same
#name rather than group by looking at the refraction index - the latter
#is more appropriate, eg., in lensParam, where I am computing physical
#parameters for each chunk of glass.

proc lensSeidel {optic xmm ymm icolor {icolor2 0}} {
   seidel $optic $xmm $ymm $icolor $icolor2
   set format "%8s %9s %9s %9s %9s"
   set format2 " %9s %9s"
   puts -nonewline stdout [format $format Lens TAS TSC TSA PETZ]
   if {$icolor2 != 0} {puts -nonewline stdout [format $format2 TAC TLC]}
   puts stdout " (microns)"
   set format "%8s %9.0f %9.0f %9.0f %9.0f"
   set format2 " %9.0f %9.0f"
   set tassum 0.
   set tscsum 0.
   set tsasum 0.
   set tacsum 0.
   set tlcsum 0.
   set hsum 0.

   set r [expr sqrt($xmm*$xmm + $ymm*$ymm)]
   set lastname "bogusname"
   set taslsum 0.
   set tsclsum 0.
   set tsalsum 0.
   set taclsum 0.
   set tlclsum 0.
   set hlsum 0.

#End limit includes focal plane, which will not get printed (unless I have
#2 surfaces with the name "FOCAL", e.g.)
   loop i 1 [exprGet $optic.diagram->np] {
	set isurf [exprGet $optic.diagram->indx<$i>]
	set surfid [surfId $optic $isurf]
	set name [showName $optic $surfid]
	set tas [expr [exprGet $optic.optic<$isurf>->tas]*1.e3]
	set tsc [expr [exprGet $optic.optic<$isurf>->tsc]*1.e3]
	set tsa [expr [exprGet $optic.optic<$isurf>->tsa]*1.e3]
	set petz [exprGet $optic.optic<$isurf>->petz]
	set tac [expr [exprGet $optic.optic<$isurf>->tac]*1.e3]
	set tlc [expr [exprGet $optic.optic<$isurf>->tlc]*1.e3]

#Express field curvature as height of focal plane at edge.
	set height [expr .5*$petz*$r*$r*1.e3]


	set tassum [expr $tassum + $tas]
	set tscsum [expr $tscsum + $tsc]
	set tsasum [expr $tsasum + $tsa]
	set tacsum [expr $tacsum + $tac]
	set tlcsum [expr $tlcsum + $tlc]
	set hsum [expr $hsum + $height]


#Typical units are microns.
	if {$name != "$lastname"} {
	   if {$i != 1} {
		puts -nonewline stdout [format $format $lastname $taslsum \
		   $tsclsum $tsalsum $hlsum]
		if {$icolor2 != 0} {
		   puts -nonewline stdout [format $format2 $taclsum $tlclsum]
		   }
		puts stdout ""
		}
	   set taslsum $tas
	   set tsclsum $tsc
	   set tsalsum $tsa
	   set taclsum $tac
	   set tlclsum $tlc
	   set hlsum $height
	   set lastname $name
	} else {
	   set taslsum [expr $taslsum + $tas]
	   set tsclsum [expr $tsclsum + $tsc]
	   set tsalsum [expr $tsalsum + $tsa]
	   set taclsum [expr $taclsum + $tac]
	   set tlclsum [expr $tlclsum + $tlc]
	   set hlsum [expr $hlsum + $height]
	   }
	}
   puts -nonewline stdout [format $format Sum $tassum $tscsum $tsasum $hsum]
   if {$icolor2 != 0} {
	puts -nonewline stdout [format $format2 $tacsum $tlcsum]
	}
   puts stdout ""
   return
   }

###################################################################
#Compute aberrations from actual ray traces.
#For now, just do lateral chromatic.  Axial is all tied up with other
#aberrations (including spherical).

proc aberrRay {hndl xmm ymm icolor1 icolor2} {

#Compute the average scale factor.
   set sum 0.
   set n 0
   foreach ifil "$icolor1 $icolor2" {
	ray $hndl $xmm $ymm 0 0 $ifil 0
	set np [exprGet $hndl.diagram->np]
	set np1 [expr $np-1]
	set x($ifil) [exprGet $hndl.diagram->xray<$np1>]
	set y($ifil) [exprGet $hndl.diagram->yray<$np1>]
	}
   set lateral [expr sqrt(pow($x($icolor1)-$x($icolor2),2) + \
	pow($y($icolor1)-$y($icolor2),2))]
   return [format %.3f $lateral]
   }

#########################################################################
#Plot aberration vs. field radius.

proc seidelPlot {optic aberr surfid ifil {ifil2 ""}} {
   set aberr [string tolower $aberr]
   set xrad [showFocal $optic $ifil xrad]

#Just go out in xrad
   set xrad [expr abs($xrad)]
   set min 1.e10
   set max -1.e10
   foreach afract "0. .2 .4 .6 .8 1." {
	set fract $afract
	set xmm [expr $xrad*$fract]
	set ymm 0.
	seidel $optic $xmm $ymm $ifil $ifil2

	loop i 1 [exprGet $optic.diagram->np] {
	   set isurf [exprGet $optic.diagram->indx<$i>]
	   set id [surfId $optic $isurf]
	   if {$id != $surfid} continue
	   set name [showName $optic $surfid]
	   set val($xmm) [expr [exprGet $optic.optic<$isurf>->$aberr]*1.e3]
	   if {$aberr == "petz"} {

#Express field curvature as height of focal plane at edge.
		set val($xmm) [expr .5*$val($xmm)*$r*$r*1.e3]
		}
	   set min [expr min($min,$val($xmm))]
	   set max [expr max($max,$val($xmm))]
	   }

#Repeat for negative xmm.
	set fract [expr -1.*$afract]
	set xmm [expr $xrad*$fract]
	set ymm 0.
	seidel $optic $xmm $ymm $ifil $ifil2

	loop i 1 [exprGet $optic.diagram->np] {
	   set isurf [exprGet $optic.diagram->indx<$i>]
	   set id [surfId $optic $isurf]
	   if {$id != $surfid} continue
	   set name [showName $optic $surfid]
	   set val($xmm) [expr [exprGet $optic.optic<$isurf>->$aberr]*1.e3]
	   if {$aberr == "petz"} {

#Express field curvature as height of focal plane at edge.
		set val($xmm) [expr .5*$val($xmm)*$r*$r*1.e3]
		}
	   set min [expr min($min,$val($xmm))]
	   set max [expr max($max,$val($xmm))]
	   }
	}
   if {$min == $max} {
	set min [expr -1.*$max]
	set max [expr 1.5*$max]
	}
   pgEnv -$xrad $xrad $min $max 0 0
   pgLabel XMM "Aberration (Microns)" [string toupper $aberr]
   set xlist ""
   set ylist ""
   foreach xmm [lsort -real [array names val]] {
	lappend xlist $xmm
	lappend ylist $val($xmm)
	}
   pgLine $xlist $ylist
   return
   }

#########################################################################
#Plot aberrations from zernikes vs. field radius

proc zradPlot {optic ifil} {
   set xrad [showFocal $optic $ifil xrad]

#Just go out in xrad
   set xrad [expr abs($xrad)]
   set min 1.e10
   set max -1.e10
   set zern [genericNew ZERNIKE]
   set wave [showWave $optic $ifil]
   set diam [telDiam $optic $ifil]
   set exit [showFocal $optic $ifil exit]
   set xmag [showFocal $optic $ifil xmag]

#Factor to convert from wave-front error to focal plane tranverse error.
#Note that I still need to multiply by a dimensionless coefficient that
#depends on each aberration type.
   set fact [expr abs($exit/($diam*$xmag))]

   foreach afract "0. .1 .2 .3 .4 .5 .6 .7 .8 .9 1." {
	set fract $afract
	set xmm [expr $xrad*$fract]
	set ymm 0.
	zernikeFit $optic $zern 4 $xmm $ymm $ifil

#Astigmatism
#	set n 2
#	set m 2
#	set N [expr $n + $m]
#	set k [expr ($n - $m)/2]
#	set j [expr ($N/2)*($N/2) + 2*$k]
#	set j 4
	set astig($xmm) [expr [exprGet $zern.coeff<4>]*$wave * 2. * $fact]
	set coma($xmm) [expr [exprGet $zern.coeff<6>]*$wave * 3. * $fact]
	set sa($xmm) [expr [exprGet $zern.coeff<8>]*$wave * 4. * $fact]
	set min [min $min $astig($xmm) $coma($xmm) $sa($xmm)]
	set max [max $max $astig($xmm) $coma($xmm) $sa($xmm)]

#Repeat for negative xmm.
	set fract [expr -1.*$afract]
	set xmm [expr $xrad*$fract]
	set ymm 0.
	zernikeFit $optic $zern 4 $xmm $ymm $ifil
	set astig($xmm) [expr [exprGet $zern.coeff<4>]*$wave * 2. * $fact]
	set coma($xmm) [expr [exprGet $zern.coeff<6>]*$wave * 3. * $fact]
	set sa($xmm) [expr [exprGet $zern.coeff<8>]*$wave * 4. * $fact]
	set min [min $min $astig($xmm) $coma($xmm) $sa($xmm)]
	set max [max $max $astig($xmm) $coma($xmm) $sa($xmm)]
	}
   zernikeDel $zern
   if {$min == $max} {
	set min [expr -1.*$max]
	set max [expr 1.5*$max]
	}
   pgSci 1
   pgSch 1
   pgEnv -$xrad $xrad $min $max 0 0
   pgLabel XMM "Transverse Aberr (Microns)" "Coma, Astigmatism, Spherical"

#Coma
   set xlist ""
   set ylist ""
   foreach xmm [lsort -real [array names coma]] {
	lappend xlist $xmm
	lappend ylist $coma($xmm)
	}
   pgSci 2
   pgLine $xlist $ylist
   pgText [expr -.2*$xrad] [expr $min + .95*($max-$min)] COMA

#Astigmatism
   set xlist ""
   set ylist ""
   foreach xmm [lsort -real [array names astig]] {
	lappend xlist $xmm
	lappend ylist $astig($xmm)
	}
   pgSci 3
   pgLine $xlist $ylist
   pgText [expr -.2*$xrad] [expr $min + .9*($max-$min)] ASTIG

#Spherical Aberration
   set xlist ""
   set ylist ""
   foreach xmm [lsort -real [array names sa]] {
	lappend xlist $xmm
	lappend ylist $sa($xmm)
	}
   pgSci 4
   pgLine $xlist $ylist
   pgText [expr -.2*$xrad] [expr $min + .85*($max-$min)] SPHERICAL

   pgSci 1
   return
   }

