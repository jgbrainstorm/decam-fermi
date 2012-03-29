# Steve Kent and al. additions/modifications to dervishStartup.tcl
#
# handleShow		Show the value of what a handle points to.
##############################################################################


proc handleShow {hndl} {
    set test [schemaKindGet $hndl]
    if {$test == "PRIM"} {

	# Handle is a primitive type.

        set type [lindex [schemaGet $hndl] 0]

	# lreplace function gets rid of the leading type returned by membersGet
	#	set val [lindex [lreplace [membersGet $hndl] 0 0] 0]
	# the new (1/95) exprGet by default does not return an header so

        set val [lindex [exprGet $hndl] 0]
        set val [string trim $val { }]

	if {$type != "char"} {
	    return $val
	}
	# character.  Return just 1 character.  Substitute non-printing charactesr
	# with a .
	if {[string length $val] == 3} then {
	    set val [string range $val 1 1]
	    return $val
	} else {
	    set val ""
	    return $val
	}
    }
    if {$test == "ENUM"} {
	
	#handle is an enumerated type

        return [exprGet -enum $hndl]
    }

    # Handle is a complicated or unknown type.

    return $hndl
}
