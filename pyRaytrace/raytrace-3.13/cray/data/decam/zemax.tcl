#Read zemax files and tweak to match my conventions.
#number is 2602 or 2605

proc rbRead {design} {
   set hndl [zemaxRead $design.zmx]
   set number 2
   opticRemove $hndl 3
   opticRemove $hndl 3
   if {$number == 3} {
	opticRemove $hndl 3
	} 
   setName $hndl 1 STOP
   setName $hndl 3 C1
   setName $hndl 4 C1
   setName $hndl 5 C2
   setName $hndl 6 C2
   setName $hndl 7 C3
   setName $hndl 8 C3
   setName $hndl 9 FILTER
   setName $hndl 10 FILTER
   setName $hndl 11 C4
   setName $hndl 12 C4
   setName $hndl 13 C5
   setName $hndl 14 C5
   setName $hndl 15 FOCAL

#zemaxRead now defines only 4 wavelengths.  Add 4 more
   waveSwitch $hndl 1 .4
   waveSwitch $hndl 2 .54
   waveSwitch $hndl 3 .56
   waveSwitch $hndl 4 .68
   waveAdd $hndl .69 1
   waveAdd $hndl .82 1
   waveAdd $hndl .82 1
   waveAdd $hndl 1.08 1

#Insert 3 surfaces for new primaries.
   set z2 [showSurf $hndl 2 z]
   set i2 [showIndex $hndl 2 1]
   set glass2 [showGlass $hndl 2]
   set stop2 [showSurf $hndl 2 outstop]
   set curv2 [showSurf $hndl 2 curv]
   set ccon2 [showSurf $hndl 2 ccon]

   setSurfId $hndl 2 2.01
   foreach i "2 3 4" {
	opticInsert $hndl 2.01 2
	}

   foreach i "1 2 3 4" {
	setName $hndl 2.0$i PRIMARY
	setGlass $hndl 2.0$i $glass2
	setSurf $hndl 2.0$i z $z2
	setSurf $hndl 2.0$i stoptype 2
	setSurf $hndl 2.0$i outstop $stop2
	setSurf $hndl 2.0$i curv $curv2
	setSurf $hndl 2.0$i ccon $ccon2
	foreach j "1 2 3 4 5 6 7 8" {
	   setIndex $hndl 2.0$i $j 0
	   }
	setIndex $hndl 2.0$i [expr 2*$i-1] $i2
	setIndex $hndl 2.0$i [expr 2*$i] $i2
	}

   foreach i "3 4 5 6 7 8 9 10 11 12 13 14 15" {
	setSurfId $hndl [expr $i+3] $i
	}

#Reset optical indexes
   foreach i "1 2 3 4 5 6 7 8" {
	set wave [showFocal $hndl $i wave]
	waveSwitch $hndl $i $wave
	}

#Inner stop
   setSurf $hndl 1 instop 700.
   setSurf $hndl 1 stoptype 1
   opticInfo $hndl
   scaleGuess $hndl
   stopcomp $hndl
   opticinc $hndl 1
   return $hndl
   }
