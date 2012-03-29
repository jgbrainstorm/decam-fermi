#Create an optical design from reading gladder's designs

#decamSpot
#opticGhost $hndl "3 4 5 6 7 8 9 10 11 12 13 14" 2 64
#ghostAdd (in kentools)
#tolExpand
#tol1Flags
#asphTol1
# ... 2 3 4
#tolRMS		compute rms in tolerance Table.
#flags	Least squares flags for tweaking.
#imgSummary
#may07
#may11
#tolCarlo <hndl> <niter>
#tempTol	Tolerance temperature
#adcTest <hndl> filter zenith	Compute image width for given SDSS filter, z
#######################################################################
#Define a spot pattern where the last spot has lscale = 1

proc lsstSpot {} {
   clearSpot
   addSpot 0 0 .17
   addSpot .2 0 .2
   addSpot .4 0 .4
   addSpot .6 0 .6
   addSpot .8 0 .8
   addSpot 1.0 0 2. 1
   return
   }

#####################################################################
#Seppala LSST design.

proc lsst {} {
   set fid [open lsst.txt]
   set line [gets $fid]
   set optic [opticNew]

#Filter wavelengths
   set waveList(1) [expr .3577-.0323]
   set waveList(2) [expr .3577+.0323]
   set waveList(3) [expr .436 - .0495]
   set waveList(4) [expr .436 + .0495]
   set waveList(5) [expr .537 - .047]
   set waveList(6) [expr .537 + .047]
   set waveList(7) [expr .644 - .0755]
   set waveList(8) [expr .644 + .0755]
   set waveList(9) [expr .8075 - .075]
   set waveList(10) [expr .8075 + .075]
   set waveList(11) [expr .940 - .1]
   set waveList(12) [expr .940 + .1]


#Focal plane
   foreach ifil [array names waveList] {
	setFocal $optic $ifil xpos 0
	setFocal $optic $ifil ypos 0
	setFocal $optic $ifil xrad 275.
	setFocal $optic $ifil yrad 275.
	setFocal $optic $ifil wave $waveList($ifil)
	setFocal $optic $ifil weight 1
	setIndex $optic 0 $ifil 1
	}
   handleSet $optic.ncolor [array size waveList]
   setSurf $optic 0 z -1.e14

   set z 0

#I will not insert stop at the top.  This way, I can add conics, etc.
#below in a way that matches the surface number in the file lsst.txt.
#I insert stop at end.
   while {1} {
	set line [gets $fid]
echo $line
	if {[eof $fid]} break
	set surf [lindex $line 0]
	set radius [lindex $line 1]
	if {$radius != 0.} {set curv [expr 1./$radius]} else {set curv 0.}
	set thick [lindex $line 2]
	set glass [lindex $line 3]
	set name [lindex $line 4]

#Primary - I will create 6 surfaces, one for each filter, and allow refocus
#for each filter.
	if {$surf == 1} {
	   foreach filter "1 2 3 4 5 6" {
		set surfid $surf.[format %02d $filter]
		setSurf $optic $surfid z $z
		setSurf $optic $surfid curv $curv

#Parameters not in file
		setSurf $optic $surfid stoptype 2
		setSurf $optic $surfid outstop 4250.
		setName $optic $surfid $name
		setGlass $optic $surfid $glass

#Hyperboloid on primary
		setSurf $optic $surfid ccon -.960331
		setSurf $optic $surfid a6 [expr -1.0354e-10*1.e-15]
		set ifil1 [expr $filter*2-1]
		set ifil2 [expr $ifil1+1]
		foreach ifil "$ifil1 $ifil2" {
		   set index [glass $glass $waveList($ifil)]
		   set index [expr -1.*$index]
		   setIndex $optic $surfid $ifil $index
		   }
		}
	} else {
	   setSurf $optic $surf z $z
	   setSurf $optic $surf curv $curv
	   foreach ifil [array names waveList] {
		set index [glass $glass $waveList($ifil)]
		set index [expr -1.*$index]
		setIndex $optic $surf $ifil $index
		}
	   setName $optic $surf $name
	   setGlass $optic $surf $glass
	   }
	set z [expr $z + $thick]
	}
   close $fid
   incr surf

#Focal plane
   setSurf $optic $surf z $z
   foreach ifil [array names waveList] {
	setIndex $optic $surf $ifil -1
	}
   setGlass $optic $surf air
   setName $optic $surf FOCAL

#Add higher order terms.
   setSurf $optic 2 ccon -0.160446
   setSurf $optic 3 ccon 0.016852
   setSurf $optic 5 ccon 2.052976
   setSurf $optic 7 ccon -0.646026

   setSurf $optic 2 a6 [expr -1.3203e-5*1.e-15]
   setSurf $optic 3 a6 [expr 4.2235e-7*1.e-15]
   setSurf $optic 5 a6 [expr 0.019629*1.e-15]
   setSurf $optic 7 a6 [expr -0.059479*1.e-15]

   setSurf $optic 2 a8 [expr 2.4052e-8*1.e-21]
   setSurf $optic 3 a8 [expr 6.5001e-9*1.e-21]
   setSurf $optic 5 a8 [expr 1.4461e-5*1.e-21]
   setSurf $optic 7 a8 [expr -0.070993*1.e-21]

   setSurf $optic 2 a10 [expr -9.7751e-8*1.e-27]

#Insert stop
   stopAdd $optic 0 -5625. 1600. 6000
   rayPattern $optic 6 1
   colorcount $optic

#Design name
   setDesign $optic lsst

#Scale factor
   set scale 0
   set n 0
   foreach ifil [array names waveList] {
	opticInfo $optic $ifil
	set fl [showFocal $optic $ifil fl]
	set scale [expr $scale + 206265./$fl]
	incr n
	}
   set scale [expr -1.*$scale/$n]

#Hardwire - paraxial calculation is not very good.
   set scale 19.5951
   foreach ifil [array names waveList] {
	setFocal $optic $ifil scale $scale
	}
   stopcomp $optic
   opticinc $optic 1
   return $optic
   }



