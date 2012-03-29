#Create an optical design from reading pfadc.txt
#Commands:
#   blanco
#   mayall

#This is the Blanco 4 m corrector
proc blanco {} {
   set fid [open pfadc.txt]
   set line [gets $fid]
   set optic [opticNew]

#Filter wavelengths
   set waveList(1) .39 
   set waveList(2) .5
   set waveList(3) .65
   set waveList(4) .90

#Focal plane
   foreach ifil [array names waveList] {
	setFocal $optic $ifil xpos 0
	setFocal $optic $ifil ypos 0
	setFocal $optic $ifil xrad 85.
	setFocal $optic $ifil yrad 85.
	setFocal $optic $ifil wave $waveList($ifil)
	setFocal $optic $ifil weight 1
	setIndex $optic 0 $ifil 1
	}
   handleSet $optic.name "CTIO-4M-PFADC"
   handleSet $optic.ncolor [array size waveList]
   setSurf $optic 0 z -1.e14

   set z 0
   while {1} {
	set line [gets $fid]
	if {[eof $fid]} break
	set surf [lindex $line 0]
	set rad [lindex $line 1]
	if {$rad != 0.} {
	   set curv [expr -1./$rad]
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
	   set index [expr -1.*$index]
	   setIndex $optic $surf $ifil $index
	   }
	setGlass $optic $surf $glass
	setName $optic $surf $name
	set z [expr $z - $thick]
	}
   incr surf
   setSurf $optic $surf z $z
   foreach ifil [array names waveList] {
	setIndex $optic $surf $ifil -1
	}
   setGlass $optic $surf air
   setName $optic $surf FOCAL

#Parameters not in file
   setSurf $optic 1 stoptype 2
   setSurf $optic 1 outstop 2000.

   setSurf $optic 1 ccon -1.09763

#Surface 3 is tilted 1.17 deg
#Surface 6 is tilted 1.37 deg
   setSurf $optic 3 theta [expr 1.17/57.3]
   setSurf $optic 6 theta [expr 1.37/57.3]

#This will increment most surface numbers by 1.
   stopAdd $optic 0 -7800 700 4000
   setGlass $optic 1 air
   setName $optic 1 STOP

   rayPattern $optic 6 1
   colorcount $optic

#Scale factor
   set sum 0
   set n 0
   foreach ifil [array names waveList] {
	opticInfo $optic $ifil
	set fl [showFocal $optic $ifil fl]
	set scale [expr 206265./$fl]
	set sum [expr $sum + $scale]
	incr n
	}
   set sum [expr 1.*$sum/$n]
   foreach ifil [array names waveList] {
	setFocal $optic $ifil scale $sum
	}
   stopcomp $optic
   opticinc $optic 1
   return $optic
   }

###############################################################33
#This is the Mayall 4 m PF corrector

proc mayall {} {
   set fid [open mayall.txt]
   set line [gets $fid]
   set optic [opticNew]

#Filter wavelengths
   set waveList(1) .39 
   set waveList(2) .5
   set waveList(3) .65
   set waveList(4) .90

#Focal plane
   foreach ifil [array names waveList] {
	setFocal $optic $ifil xpos 0
	setFocal $optic $ifil ypos 0
	setFocal $optic $ifil xrad 85.
	setFocal $optic $ifil yrad 85.
	setFocal $optic $ifil wave $waveList($ifil)
	setFocal $optic $ifil weight 1
	setIndex $optic 0 $ifil 1
	}
   handleSet $optic.name "MAYALL-4M-PFADC"
   handleSet $optic.ncolor [array size waveList]
   setSurf $optic 0 z -1.e14

   set z 0
   while {1} {
	set line [gets $fid]
	if {[eof $fid]} break
	set surf [lindex $line 0]
	set rad [lindex $line 1]
	if {$rad != 0.} {
	   set curv [expr -1./$rad]
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
	   set index [expr -1.*$index]
	   setIndex $optic $surf $ifil $index
	   }
	setGlass $optic $surf $glass
	setName $optic $surf $name
	set z [expr $z - $thick]
	}
   incr surf
   setSurf $optic $surf z $z
   foreach ifil [array names waveList] {
	setIndex $optic $surf $ifil -1
	}
   setGlass $optic $surf air
   setName $optic $surf FOCAL

#Parameters not in file
   setSurf $optic 1 stoptype 2
   setSurf $optic 1 outstop 2000.
   setSurf $optic 1 ccon -1.09763

#Angles on Risley prism surfaces (in order from primary to focal plane):
# .076 deg
# 2.100 deg
# 0
# -0.074 deg
# -2.233
# 0
#The non-zero tilts on the outer surfaces are to compensate for the fact that
#UBK7 and LLF6 are not perfectly matched at the median wavelength.

#Hmm, this is not quite right and does not work yet.
   setSurf $optic 8 theta [expr .076/57.3]
   setSurf $optic 9 theta [expr 2.1/57.3]
   setSurf $optic 11 theta [expr -.074/57.3]
   setSurf $optic 12 theta [expr -2.233/57.3]

#This will increment most surface numbers by 1.
   stopAdd $optic 0 -7800 700 4000
   setGlass $optic 1 air
   setName $optic 1 STOP

   rayPattern $optic 6 1
   colorcount $optic

#Scale factor
   set sum 0
   set n 0
   foreach ifil [array names waveList] {
	opticInfo $optic $ifil
	set fl [showFocal $optic $ifil fl]
	set scale [expr 206265./$fl]
	set sum [expr $sum + $scale]
	incr n
	}
   set sum [expr 1.*$sum/$n]
   foreach ifil [array names waveList] {
	setFocal $optic $ifil scale $sum
	}
   stopcomp $optic
   opticinc $optic 1
   return $optic
   }

##########################################################################
#Input radii of primary, secondary mirror and separation.  Compute all 
#other RC parameters.  This is easier than using secondary conic in place 
#of R1 and trying to go backwards (messy eqns.)
#
#The goal is (was) to see if one could infer the proper primary mirror
#parameters from the known secondary mirror params.  However, we don't get
#total consistency assuming that the F/7.8 mirror is a true RC secondary.
#Also, the other secondaries couldn't be RC secondaries unless the 
#backfocal distance is something totally unreasonable.


proc rc {r1 r2 sep} {

   set f1 [expr -$r1/2.]
   set f2 [expr -$r2/2.]
   set rho [expr $r2/$r1]
   echo rho $rho
   set k [expr 1.-$sep/$f1]
   echo k $k
   set m [expr $rho/($rho-$k)]
   echo mag $m
   set beta [expr $k*($m+1.)-1.]
   echo beta $beta

#Predicted conics
   set k1 [expr -1.-2.*(1.+$beta)/($m*$m*($m-$beta))]
   set k2 [expr -pow(($m+1.)/($m-1.),2) - \
	2.*$m*($m+1.)/(($m-$beta)*pow($m-1,3))]
   echo conic1: Predicted [format %.4f $k1]
   echo conic2: Predicted [format %.4f $k2]
   return
   }
