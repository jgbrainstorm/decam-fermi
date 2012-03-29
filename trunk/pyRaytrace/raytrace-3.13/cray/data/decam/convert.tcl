#Convert blanco-2602-v203 to a full design with different filter locations
#and more accurate stop

proc convert {} {
   set hndl [lensRead blanco-2602-v203]

   echo Resetting filter surf ids
   setSurfId $hndl 9 9.01
   setSurfId $hndl 10 10.01
   echo Adding new filters
   surfInsert $hndl 9.01 [showSurf $hndl 9.01 z] FSEPPS 2
   surfInsert $hndl 9.02 [showSurf $hndl 9.01 z] FSEPPS 2
   surfInsert $hndl 9.03 [showSurf $hndl 9.01 z] FSEPPS 2
   surfInsert $hndl 10.01 [showSurf $hndl 10.01 z] air 2
   surfInsert $hndl 10.02 [showSurf $hndl 10.01 z] air 2
   surfInsert $hndl 10.03 [showSurf $hndl 10.01 z] air 2

#Adjust filter positions
   echo Adjusting filter heights
   incrSurf $hndl 9.02 z 38
   incrSurf $hndl 9.03 z 76
   incrSurf $hndl 9.04 z 114
   incrSurf $hndl 10.02 z 38
   incrSurf $hndl 10.03 z 76
   incrSurf $hndl 10.04 z 114
   echo Adjusting filter colors
   keepIndex $hndl 9.01 1 2
   keepIndex $hndl 9.02 3 4
   keepIndex $hndl 9.03 5 6
   keepIndex $hndl 9.04 7 8
   keepIndex $hndl 10.01 1 2
   keepIndex $hndl 10.02 3 4
   keepIndex $hndl 10.03 5 6
   keepIndex $hndl 10.04 7 8

#Create a stop
   echo New stops
   set rad 1967
   set index [surfIndex $hndl 2.01]
   set zoff [zsurf $hndl $index $rad 0]

   surfInsert $hndl 1 $zoff air 0 APSTOP

#Reset stop types of primary mirror
   foreach i "1 2 3 4" {
	setSurf h0 3.0$i stoptype 0
	}
   setSurf $hndl 2 stoptype 2
   setSurf $hndl 2 outstop $rad

#Refocus
   echo Refocus
   foreach i "1 2 3 4" {
	set i1 [expr $i+1]
	fastLstsqInit
	linkColorFlag h0 "$i $i1" 2
	setSurfFlag h0 3.0$i z 1
	decamSpot
	fastLstsq $hndl
	}

#Reset scale so I hit the edge of the focal plane
   echo New scale
   clearFlags
   fastLstsqInit 1
   verbose 0
   linkColorFlag $hndl "1 2 3 4 5 6 7 8" 2
   linkFocalFlag $hndl "1 2 3 4 5 6 7 8" scale 2
   linkFocalFlag $hndl "1 2 3 4 5 6 7 8" dist 4
   decamSpot
   fastLstsq $hndl

   echo New stops
   stopComp $hndl
   opticInc $hndl
   return $hndl
   }
