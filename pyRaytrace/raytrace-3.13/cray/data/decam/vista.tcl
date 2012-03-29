#Create an optical design from subaru.txt

#This is the VISTA 4 m optical corrector.
proc vista {} {
   set fid [open vista.txt]
   set line [gets $fid]
   set optic [opticNew]

#Filter wavelengths
   set waveList(1) .39 
   set waveList(2) .54
   set waveList(3) .56
   set waveList(4) .68
   set waveList(5) .69
   set waveList(6) .82
   set waveList(7) .82
   set waveList(8) 1.08

#Focal plane
   foreach ifil [array names waveList] {
	setFocal $optic $ifil xpos 0
	setFocal $optic $ifil ypos 0
	setFocal $optic $ifil xrad 100.
	setFocal $optic $ifil yrad 100.
	setFocal $optic $ifil xrad 225.3
	setFocal $optic $ifil yrad 225.3
	setFocal $optic $ifil wave $waveList($ifil)
	setFocal $optic $ifil weight 1
	setIndex $optic 0 $ifil 1
	}
   handleSet $optic.name "VISTA OPTICAL CORRECTOR"
   handleSet $optic.ncolor [array size waveList]
   setSurf $optic 0 z -1.e14

   set z 0
   while {1} {
	set line [string trim [gets $fid]]
	if {[eof $fid]} break
	if {[string length $line] == 0} break
	set surf [lindex $line 0]
	set rad [lindex $line 1]
	if {$rad != 0.} {
	   set curv [expr 1./$rad]
	} else {
	   set curv 0.
	   }
	set thick [lindex $line 2]
	set glass [lindex $line 3]
	set name [lindex $line 4]
	setSurf $optic $surf z $z
	setSurf $optic $surf curv $curv
	foreach ifil [array names waveList] {
	   set index [glass $glass $waveList($ifil)]
	   setIndex $optic $surf $ifil $index
	   }
	setGlass $optic $surf $glass
	setName $optic $surf $name
	set z [expr $z + $thick]
	}
   incr surf
   setSurf $optic $surf z $z
   setName $optic $surf FOCAL
   foreach ifil [array names waveList] {
	setIndex $optic $surf $ifil 1
	}

#Parameters not in file
   setSurf $optic 1 stoptype 2
   setSurf $optic 1 outstop 2005.
   setSurf $optic 1 ccon -1.129792

   setSurf $optic 2 ccon -5.548792

#Tilts
   setSurf $optic 6 theta [expr -.56/57.3]
   setSurf $optic 7 theta [expr -.018/57.3]

#I reverse the sign of the next angles from the default design
#rather than set phi = pi.
   setSurf $optic 9 theta [expr .59/57.3]
   setSurf $optic 10 theta [expr .019/57.3]

   setSurf $optic 14 a4 -9.6638541e-11
   setSurf $optic 14 a6 -3.6099989e-15
   setSurf $optic 14 a8 -1.6380462e-22
   setSurf $optic 14 a10 -3.6121901e-25

   stopAdd $optic 0 -2725. 620.3 2200.
   setName $optic 1 STOP
   setGlass $optic 1 air
   rayPattern $optic 6 1
   colorcount $optic

#Scale factor.  Use a single value for all filters.
   set sum 0.
   set n 0
   foreach ifil [array names waveList] {
	opticInfo $optic $ifil
	set fl [showFocal $optic $ifil fl]
	set scale [expr 206265./$fl]
	set sum [expr $sum + $scale]
	incr n
	}
   set scale [expr $sum/$n]

#Better yet (since there is distortion)
   set scale 16.90556
   foreach ifil [array names waveList] {
	setFocal $optic $ifil scale $scale
	}
   stopcomp $optic
   opticinc $optic 1

#Secondary is undersized

#I can also undersize the primary
#This is much nicer, because I can use a coarser rayPattern and get
#accurate image sizes.  With a properly undersized secondary, I need to
#use rayPattern 20, which is SSLLOOWW
   setSurf $optic 2 outstop 1850.

   stopcomp $optic
#   setSurf $optic 3 stoptype 1
#   setSurf $optic 3 outstop 620.3

   return $optic
   }



