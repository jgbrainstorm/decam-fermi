#Analyze spectrograph design

#Flags and spots.  I will designate 1 spot to be a constraint on scale
#To constrain at only one wavelength, I will allow scale to be free on
#all other wavelengths.

#Scale factor is basically 1./camera-focal-length.

#Flags for untilted configuration.
proc specFlags {hndl} {
   clearFlags
   clearSpot
   addSpot 0 0 1 0
   addSpot .7 0 1 0
   addSpot 1 0 1 1
   addSpot .7 .7 1 0
   addSpot 1 .7 1 0
   addSpot 1 1 1 0
   set ncolor [exprGet $hndl.ncolor]
   set ifils ""
   foreach ifil [range 1-$ncolor] {
	lappend ifils $ifil
	}

#Scale factor is negative
#Set xrad and yrad negative
   foreach ifil $ifils {
	set scale [showFocal $hndl $ifil scale]
	set scale [expr -1.*abs($scale)]
	setFocal $hndl $ifil scale $scale
	set xrad [showFocal $hndl $ifil xrad]
	set xrad [expr -1.*abs($xrad)]
	setFocal $hndl $ifil xrad $xrad
	set yrad [showFocal $hndl $ifil yrad]
	set yrad [expr -1.*abs($yrad)]
	setFocal $hndl $ifil yrad $yrad
	}
   foreach ifil $ifils {
	setColorFlag $hndl $ifil 1
	}

#Keep scale factor fixed for filter 3.  This wavelength falls near the
#middle of the spectrum.
   foreach ifil $ifils {
	if {$ifil == 2} continue
	setFocalFlag $hndl $ifil scale 1
	}

   setSurf $hndl 17 a4 0
   setSurf $hndl 17 a6 0
   setSurf $hndl 17 a8 0
#   setSurfFlag $hndl 17 a4 1
#   setSurfFlag $hndl 17 a6 1
#   setSurfFlag $hndl 17 a8 1

   setSurfFlag $hndl 6 curv 1
   setSurfFlag $hndl 7 curv 1

   setSurfFlag $hndl 8 curv 1
   linkSurfFlag $hndl "9 10" curv 10
   linkSurfFlag $hndl "11 12" curv 12
   setSurfFlag $hndl 13 curv 1

   setSurfFlag $hndl 14 curv 1
#   linkSurfFlag $hndl "14 15 16" curv 14
   setSurfFlag $hndl 15 curv 1
   linkSurfFlag $hndl "16 17" curv 16
#   setSurfFlag $hndl 16 curv 1
#   setSurfFlag $hndl 17 curv 1
   linkSurfFlag $hndl "18 19" curv 18
#   setSurfFlag $hndl 18 curv 1
#   setSurfFlag $hndl 19 curv 1

   for {set i 6} {$i <= 20} {incr i} {
	setSurfInc $hndl $i z 1
	setSurfInc $hndl $i curv 2.e-7
	}

#   linkSurfFlag $hndl "8 9 10 11 12 13" z 2
   linkSurfFlag $hndl "14 15" z 4
   linkSurfFlag $hndl "16 17 18 19 20" z 6

   return
   }

##################################################################
#Flags to solve for scale factor only

proc scaleFlags {hndl} {

#One spot only.
   clearSpot
   addSpot 0 1. 1 1
   set ncolor [exprGet $hndl.ncolor]
   set ifils ""
   loop i 0 $ncolor {
	lappend ifils [expr $i+1]
	}

#Scale factor is negative
#Set xrad and yrad negative
   foreach ifil $ifils {
	setFocal $hndl $ifil map 2
	set scale [showFocal $hndl $ifil scale]
	set scale [expr -1.*abs($scale)]
	setFocal $hndl $ifil scale $scale
	set xrad [showFocal $hndl $ifil xrad]
	set xrad [expr -1.*abs($xrad)]
	setFocal $hndl $ifil xrad $xrad
	set yrad [showFocal $hndl $ifil yrad]
	set yrad [expr -1.*abs($yrad)]
	setFocal $hndl $ifil yrad $yrad
	}

   foreach ifil $ifils {
	setColorFlag $hndl $ifil 1
	}

#Allow xoff to float.  y direction contrains the scale factor.
#NEW: Not required in spectro mode.
#   setFocalFlag $hndl 3 xoff 1

#Keep scale factor fixed for filter 2
   foreach ifil $ifils {
	if {$ifil == 2} continue
	setFocalFlag $hndl $ifil scale 1
	}

   return
   }

##################################################################
#Flags for tilted configuration.

proc tiltFlags {hndl} {
   clearFlags
   clearSpot

   scaleFlags $hndl

#Add more spots
   addSpot 0 0 1 0
   addSpot 0 .5 1 0

#Temporary - spots for same wave solution
#   addSpot .5 0 1 0
#   addSpot -.5 0 1 0
#   addSpot .5 1 1 0
#   addSpot -.5 1 1 0

   set nsurf [exprGet $hndl.nsurf]

#Set translation increments to allow refocus at proper angle
   for {set i 6} {$i <= $nsurf} {incr i} {
	set ang [showSurf $hndl $i theta]
	setSurfInc $hndl $i z 1.
	setSurfInc $hndl $i x [expr tan($ang)]
	}

#   setSurf $hndl 15 a4 0
#   setSurf $hndl 15 a6 0
#   setSurf $hndl 15 a8 0

   setSurfFlag $hndl 15 a4 1
   setSurfFlag $hndl 15 a6 1
   setSurfFlag $hndl 15 a8 1

#   setSurfFlag $hndl 6 curv 1
#   setSurfFlag $hndl 7 curv 1
#   linkSurfFlag h0 "6 7" curv 9

   setSurfFlag $hndl 8 curv 1
   linkSurfFlag $hndl "9 10" curv 10
   linkSurfFlag $hndl "11 12" curv 12
   setSurfFlag $hndl 13 curv 1

   setSurfFlag $hndl 14 curv 1
   setSurfFlag $hndl 15 curv 1

   linkSurfFlag $hndl "16 17" curv 16
#   setSurfFlag $hndl 18 curv 1
#   setSurfFlag $hndl 19 curv 1
   linkSurfFlag $hndl "18 19" curv 18
#   setSurfFlag $hndl 18 curv 1
#   setSurfFlag $hndl 19 curv 1

#Refocus with linked x and z
   foreach i "8 9 10 11 12 13" {
#	setSurfFlag $hndl $i z 2
#	setSurfFlag $hndl $i x 2
	}

   foreach i "14 15" {
	setSurfFlag $hndl $i z 6
	setSurfFlag $hndl $i x 6
	}

   foreach i "16 17 18 19 20" {
#	setSurfFlag $hndl $i z 4
#	setSurfFlag $hndl $i x 4
	}

   foreach i "16 17 18 19 20" {
	setSurfFlag $hndl $i z 8
	setSurfFlag $hndl $i x 8
	}

#Adjust curvature of focal plane
#   setSurfFlag $hndl 20 curv 1

   return
   }

###########################################################
#Increment thickness - this means translating in both x and z to preserve
#alignment along tilted optical axis.  I actually move just 1 surface - not
#everything.

proc incrThick {hndl surfid dz} {
   set theta [showSurf $hndl $surfid theta]
   incrSurf $hndl $surfid z [expr $dz*cos($theta)]
   incrSurf $hndl $surfid x [expr $dz*sin($theta)]
   return
   }

####################################################################
#Procedure to start with spectro1u and apply tilts and set up grating
#params.  This one is a straight-through spectrograph.

proc spec1u {h} {
   lensRotate $h 3 -20.25
   incrSurf $h 3 z 14.0
   lensRotate $h 5 20.25
   incrSurf $h 5 z -14.0
   setSurf $h 4 lines 457
   setSurf $h 4 order -1
   return
   }

####################################################################
#Procedure to start with spectro2u and apply tilts and set up grating
#params.
#I believe spectro2u is designed with the following params:
#vph .55 1.05 100. 2.8 3500. 7.
#100 microns diam fiber
#R = 3500
#bragg = 7 deg

proc spec2u {h} {
#   set bragg 7.
#   set tilt [expr 57.3*asin(sin($bragg/57.3)*1.51)]
   set tilt 10.6
   lensRotate $h 3 $tilt
   surfBreak $h 3 1 $tilt
   surfBreak $h 5 1 $tilt
   setSurf $h 4 lines 457
   setSurf $h 4 order -1
   return
   }

####################################################################
#Procedure to start with spectro2u and apply tilts and set up grating
#params.
#I believe spectro2u is designed with the following params:
#vph .55 1.05 100. 2.8 3500. 7.

#spec3u is: vph .55 1.05 120 2.8 3500 6
#120 microns diam fiber
#R = 3500
#bragg = 6.5 deg

#Reduced wavelength range: vph .6 1.0 120 2.8 3500 6.5
#Reduced wavelength range: vph .6 1.0 100 2.8 3500 7.

proc spec3u {h} {
   set bragg 7.
   set tilt [expr 57.3*asin(sin($bragg/57.3)*1.51)]
   lensRotate $h 3 $tilt
   surfBreak $h 3 1 $tilt
   surfBreak $h 5 1 $tilt
   setSurf $h 4 lines 457.
   setSurf $h 4 order -1

   waveSwitch $h 1 .6
   waveSwitch $h 4 1.
   return
   }

####################################################################
#Procedure to start with spectro2ub and apply tilts and set up grating
#params.  This is blue channel of a 2-armed spectrograph.
#We keep the beam size around 158 mm.
#Camera FL is 339.1 mm
#Tilt angle:  Let "bragg" be the Bragg angle in my vph calculation.
#Let n be the refractive index of my glass.  Then sin(tilt) = n*sin(bragg)
#Fiber size is 120 microns (?) = 2.1 arcsec
#Collimator should be 439.9 mm (but my design has 427?)
#f/2.16 out

proc spec2ub {h} {
   set bragg 7.
   set n 1.514
   set tilt [expr 57.3*asin($n*$bragg/57.3)]
   set lines 562
   lensRotate $h 3 $tilt
   surfBreak $h 3 1 $tilt
   surfBreak $h 5 1 $tilt
   setSurf $h 4 lines $lines
   setSurf $h 4 order -1
   return
   }
####################################################################
#Procedure to start with spectro2ur and apply tilts and set up grating
#params.  This is red channel of a 2-armed spectrograph.
#We keep the beam size around 158 mm.
#Camera FL is 453 mm

proc spec2ur {h} {
   set bragg 7.
   set n 1.508
   set tilt [expr 57.3*asin($n*$bragg/57.3)]
   set lines 402
   lensRotate $h 3 $tilt
   surfBreak $h 3 1 $tilt
   surfBreak $h 5 1 $tilt
   setSurf $h 4 lines $lines
   setSurf $h 4 order -1
   return
   }
