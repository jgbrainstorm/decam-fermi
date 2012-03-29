# Steve Kent additions/modifications to dervishStartup.tcl
#
# chainInit		Initialise an iterator for a linked list
# chainNext		Advance an iterator
# chainPrev		Retreat an iterator
# chainList		List all handles to chains and their type
###############################################################################
#
# Chain operations
#
###############################################################################
#
# Return an iteration block as a keyed list
#
proc chainInit {chain iterator} {
	upvar 1 $iterator iter
	set cursor [chainCursorNew $chain]
	keylset iter "CHAIN" $chain "CURSOR" $cursor "HANDLE" ""
   }
#
# given an iterator block, advance it to the next element of the chain
#
proc chainNext { iterator } {
	upvar 1 $iterator iter
	set chain [keylget iter "CHAIN"]
	set cursor [keylget iter "CURSOR"]
	set handle [keylget iter "HANDLE"]
	if {$handle != ""} then {handleDel $handle}
	if {[chainSize $chain] != 0} then {
		set handle [chainWalk $chain $cursor NEXT]
	} else {
		set handle ""
		}
# Are we at end of list?
	if {$handle == ""} then {
		set iter ""
		chainCursorDel $chain $cursor
		return ""
		}
	set iter ""
	keylset iter "CHAIN" $chain "CURSOR" $cursor "HANDLE" $handle
	return $handle
}

#
# given an iterator block, advance it to the previous element of the chain
#
proc chainPrev { iterator } {
	upvar 1 $iterator iter
	set chain [keylget iter "CHAIN"]
	set cursor [keylget iter "CURSOR"]
	set handle [keylget iter "HANDLE"]
	if {$handle != ""} then {handleDel $handle}
	if {[chainSize $chain] != 0} then {
		set handle [chainWalk $chain $cursor PREVIOUS]
	} else {
		set handle ""
		}
# Are we at end of list?
	if {$handle == ""} then {
		set iter ""
		chainCursorDel $chain $cursor
		return ""
		}
	set iter ""
	keylset iter "CHAIN" $chain "CURSOR" $cursor "HANDLE" $handle
	return $handle
   }

#Remove the chain link at the current position.
proc chainPurge { iterator } {
	upvar 1 $iterator iter
	set chain [keylget iter "CHAIN"]
	set cursor [keylget iter "CURSOR"]
	set iter ""
	set handle [chainElementRemByCursor $chain $cursor]
	keylset iter "CHAIN" $chain "CURSOR" $cursor "HANDLE" ""
	return $handle
}

############################################################################
proc chainList {} {
   set list [handleListFromType CHAIN]
   foreach hndl $list {
        echo $hndl [chainTypeDefine $hndl]
        }
   }

############################################################################
proc chainAdd {chain args} {
   foreach arg $args {
	chainElementAddByPos $chain $arg
	handleDel $arg
	}
   }

############################################################################
proc chainCreate {type n} {
   set chain [chainNew $type]
   loop i 0 $n {
	chainAdd $chain [genericNew $type]
	}
   return $chain
   }

############################################################################
proc chainEach {inhndl chain fnc} {
   upvar 1 $inhndl hndl
   chainInit $chain iter
   while {1} {
	set hndl [chainNext iter]
	if {$hndl == ""} break
	eval $fnc
	}
   return
   }

#########################################################################
#Get exactly one element from chain.  Very simple predicate
proc chainGetOne {chain arg val} {
   set subchain [chainSearch $chain "{$arg = $val}"]
   if {[chainSize $subchain] == 1} {
	set hndl [chainElementGetByPos $subchain 0]
	chainDel $subchain
	return $hndl
	}
   chainDel $subchain
   return
   }

#########################################################################
#Get first matching element from chain.  Very simple predicate
proc chainGetFirst {chain arg val} {
   set subchain [chainSearch $chain "{$arg = $val}"]
   if {[chainSize $subchain] >= 1} {
	set hndl [chainElementGetByPos $subchain 0]
	chainDel $subchain
	return $hndl
	}
   chainDel $subchain
   return
   }

#########################################################################
proc chainEach {_hndl chain proc} {
   upvar 1 $_hndl hndl
   chainInit $chain iter
   while {1} {
	set hndl [chainNext iter]
	if {$hndl == ""} break
	uplevel 1 $proc
	}
   return
   }

##############################################################
proc chainGet {chain index} {
   return [chainElementGetByPos $chain $index]
}
