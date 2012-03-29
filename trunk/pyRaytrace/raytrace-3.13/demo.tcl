proc commands {} {
#Initialize display
#This is now done in rayStartup.
#plotInit a
raise .r
#Read in SNAP optical design
set hndl [opticRead snapfull]
#Draw optics layout
scaleSet 0; opticPlot $hndl 0 10 1
#Display parameters for position 11
opticDisplay $hndl 11
#Plot spot diagram for center of position 11
scaleSet .2; spotPlot $hndl 0 0 11; scaleSet 0
#Plot full focal plane layout
focalPlot $hndl 1 55
#Plot imaging array focal plane layout
focalPlot $hndl 11 55
#Plot wavefront error surface for center of position 11
waveFrontMap3d $hndl 0 0 11
#Projected view of point spread function
psfMap3d $hndl 0 0 11
#PSF's plotted at different positions
scaleSet .2; manyPlot $hndl 1; scaleSet 0
}

proc demo {} {
   set commands [split [info body commands] \n]
   lvarpop commands

   echo Demo of Ray Trace Program
   echo
   flush stdout
   after 600
   puts stdout "<press enter to proceed>: " nonewline
   flush stdout
   set char [getLine ""]
   if {$char != ""} break

   while {[llength $commands] > 0} {
	update
	echo
	set command [lvarpop commands]
	if {[string length $command] == 0} continue
	flush stdout
	after 500
	echo ---------------------------------------------------------------
	echo [string range $command 1 end]
	flush stdout
	after 1000
	set command [lvarpop commands]
	vtBold; flush stdout; echo $command; vtNormal; flush stdout
	echo
	flush stdout
	after 1000
	eval $command
	update
	echo
	flush stdout
	after 300
	puts stdout "<press enter for next step>: " nonewline
	flush stdout
	set char [getLine ""]
	if {$char != ""} break
	}
   if {[info exists hndl]} {
	echo Handle to optics data structure is $hndl
	}
   echo
   echo End of Demo
   return
   }

demo
