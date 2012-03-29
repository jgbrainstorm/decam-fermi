#Useful utilities, copied from kentools
#############################################################
#Formatted list of commands that match a pattern.
proc listcom {{pattern "*"}} {
   set comlist [info commands $pattern]
   set newcomlist [lsort $comlist]
   while {$newcomlist != ""} {
	set line ""
	for {set i 0} {$i < 3} {incr i} {
	   if {$newcomlist == ""} break
	   set command [lvarpop newcomlist]
	   append line [format "%-26s" $command]
	   }
	echo $line
	}
   }

#################################################################
#List of all TCL variables that match a pattern.
proc listvar {{pattern "*"}} {
   set comlist [info global $pattern]
   set newcomlist [lsort $comlist]
   while {$newcomlist != ""} {
	set line ""
	for {set i 0} {$i < 5} {incr i} {
	   if {$newcomlist == ""} break
	   set command [lvarpop newcomlist]
	   append line [format "%-16s" $command]
	   }
	echo $line
	}
   }

#############################################################
#Helper proc.  Write out TCL array to a file.  Allow name of file as
#an optional parameter.
   proc arrayWrite {_inarray {file ""}} {
   upvar 1 $_inarray array
   if {$file == ""} {set file $_inarray.tcl}
   set fid [open $file w]
   puts $fid "if \{\[info exists $_inarray\]\} \{unset $_inarray\}"
   foreach name [array names array] {
	puts $fid "[list set ${_inarray}($name)  $array($name)]"
	}
   puts $fid return
   close $fid
   return
   }

##############################################################
#Copy an array

proc arrayCopy {_inarray _outarray} {
   upvar 1 $_inarray inarray
   upvar 1 $_outarray outarray
   if {[info exists outarray]} {unset outarray}
   array set outarray [array get inarray]
   return
   }

#############################################################
#Assemble a message for printing.

proc msgInit {} {
   global MSG
   set MSG ""
   return
   }

#One should use msgAdd in the same way one would use the "echo" command -
#multiple arguments are added to the output line separated by spaces.
#The whole line is terminated with a <LF>

proc msgAdd {args} {
   global MSG
   foreach arg $args {
	append MSG $arg " "
	}
   append MSG \n
   return
   }

proc msgGet {} {
   global MSG
   return $MSG
   }
