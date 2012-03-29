#######################################################################
#Divide focal plane into a mosaic of detectors.
#Design originally has 4 filters.  I will add 6 dummys, then start
#mosaic numbering from 11.

proc blancoDesign {optic} {

#Reference filter
   set rfil 3

   set wave [showFocal $optic $rfil wave]
   set scale [showFocal $optic $rfil scale]

   waveAdd $optic $wave $rfil
   waveAdd $optic $wave $rfil
   waveAdd $optic $wave $rfil
   waveAdd $optic $wave $rfil
   waveAdd $optic $wave $rfil
   waveAdd $optic $wave $rfil

#Add the individual detectors in filters 11 to 18
#Number of surfaces now - must fill in optical index for new filters

   set focalid 17
   set z [showSurf $optic $focalid z]

#Optical CCDs
   set ccdwidth [expr .015*2048.]
   set ccdheight [expr .015*4096.]

#Gaps between CCDs - this is approx 114 pixels.
   set xgap 1.7
   set ygap 1.7

   set n(0) 4

#2 sides of array
   set nfil 10

#Index of focal plane
   set surfid $focalid
   foreach side "-1. 1." {

#Loop through 1 row of CCDs
	loop j 0 1 {
	   set yc [expr ($ccdheight*($j+0.5) + $j*$ygap)*$side]
	   loop i 0 $n($j) {
		set xc [expr (-($n($j)-1.)/2. + $i)*($ccdwidth + $xgap)]
		incr nfil

#Use opticInsert to add a new focal plane surface.  Get new surfid.
		set newid [opticInsert $optic $surfid]
		setName $optic $newid FOCAL
		setSurf $optic $newid x $xc
		setSurf $optic $newid y $yc
		setSurf $optic $newid z $z
		waveAdd $optic $wave $rfil
		setIndex $optic $focalid $nfil 0
		setIndex $optic $newid $nfil [showIndex $optic $focalid $rfil]
		set surfid $newid

#Regression test - ncolor should be the same as nfil
		set ncolor [exprGet $optic.ncolor]
		if {$nfil != $ncolor} {
		   error "Gak! filter count off nfil = $nfil, ncolor = $ncolor"
		   }

#Position on sky.  Don't know if my signs are correct yet!
		set xoff [expr $xc*$scale/60.]
		set yoff [expr $yc*$scale/60.]
		setFocal $optic $nfil xoff $xoff
		setFocal $optic $nfil yoff $yoff
		setFocal $optic $nfil xsize [expr -$ccdwidth/2.]
		setFocal $optic $nfil ysize [expr -$ccdheight/2.]
		setFocal $optic $nfil scale $scale
		setFocal $optic $nfil dist 0
		setFocal $optic $nfil rot 0
		opticInfo $optic $nfil
		}
	   }
	}

   stopcomp $optic
   opticinc $optic 1
   return
   }


