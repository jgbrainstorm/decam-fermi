###################################################################
#Test forking in tclX

proc forkTest {} {
   set command "Running forkEach"
   forkEach i "0 2 4 6 8 10 12 14" {
	set fid [open junk.$i w]
	puts $fid Hello
	close $fid
	puts stdout "$command $i"
	}

   set command "Running forkLoop"
   forkLoop j 0 7 {
	set fid [open junkloop.$j w]
	puts $fid Hello
	close $fid
	puts stdout "$command $j"
	}

   return
   }

##################################################################
proc forkLoop {_index begin end statement} {
   global FORK
    if {![info exists FORK]} {set FORK 1}
   upvar 1 ${_index} index
   set pidList ""
   loop index $begin $end {
	if {[info command update] != ""} {update}
	if {[llength $pidList] >= $FORK} {
	   set pid [lindex [wait] 0]
	   set lindex [lsearch $pidList $pid]
	   if {$lindex < 0} {error "Weird pid: $pid!"}
	   set pidList [lreplace $pidList $lindex $lindex]
	   }
	flush stdout
	set pid [fork]
	if {$pid != 0} {
	   lappend pidList $pid
	   continue
	} else {
	   close stdin
	   }
#We execute the statement in a catch process, which seems to catch
#any "continue" or "break" statements from executing.  This avoids our 
#having to do too much checking of unwanted behavior.
	catch {uplevel 1 $statement}
	flush stdout
	kill [pid]
	}
   foreach pid $pidList {
	wait $pid
	}
   return
   }

##################################################################
proc forkEach {_forvar forvars statement} {
   global FORK
    if {![info exists FORK]} {set FORK 1}
   upvar 1 ${_forvar} forvar
   set pidList ""
   set nproc $FORK
   foreach forvar $forvars {
	if {[info command update] != ""} {update}
	if {[llength $pidList] >= $nproc} {
	   set pid [lindex [wait] 0]
	   set lindex [lsearch $pidList $pid]
	   if {$lindex < 0} {error "Weird pid: $pid!"}
	   set pidList [lreplace $pidList $lindex $lindex]
	   }
	flush stdout
	set pid [fork]
	if {$pid != 0} {
	   lappend pidList $pid
	   continue
	} else {
	   close stdin
	   }
#We execute the statement in a catch process, which seems to catch
#any "continue" or "break" statements from executing.  This avoids our 
#having to do too much checking of unwanted behavior.
	catch {uplevel 1 $statement}
	flush stdout
	kill [pid]
	}
   foreach pid $pidList {
	wait $pid
	}
   return
   }


	