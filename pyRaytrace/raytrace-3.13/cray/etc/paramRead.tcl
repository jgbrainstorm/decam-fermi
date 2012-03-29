set tcl_precision 17

#############################################################################
# Helper proc

proc paramReadTypedef {line fid} {
   upvar 1 hdr hdr
   upvar 1 table table
   upvar 1 enum enum
   if {[llength $line] != 4} {
	close $fid
	error "Invalid typedef statement:\n$line"
	}
   set typeType [string tolower [lindex $line 1]]

#Make typeName lower case.
   set typeName [string tolower [string trim [lindex $line 3] \;]]
   if {$typeType == "enum"} {
	lappend enum() $typeName

#Trim trailing commas
	set enumList [split [string trimright [lindex $line 2] ,] ,]
	set enum($typeName) ""
	foreach val $enumList {
	   lappend enum($typeName) [string trim $val]
	   }
   } elseif {$typeType == "struct"} {
	lappend table() $typeName
	set table($typeName) ""
	set paramList [split [lindex $line 2] \;]
	foreach paramSet $paramList {
	   set paramSet [string trim $paramSet]
	   if {[llength $paramSet] < 2} continue
#Convert dimension delimiters to angle brackets. Also get rid of embedded
#whitespace.
#I could conceivably convert some enum info unintentionally, but tough.
	   regsub -all {\[[ ]*([0-9]+)[ ]*\]} $paramSet {<\1>} paramSet
	   regsub -all {\([ ]*([0-9]+)[ ]*\)} $paramSet {<\1>} paramSet
	   regsub -all {<[ ]*([0-9]+)[ ]*>} $paramSet {<\1>} paramSet

#Isolate dimension info with spaces
	   regsub -all {([^ ])<} $paramSet {\1 <} paramSet
	   regsub -all {>([^ ])} $paramSet {> \1} paramSet

#Basic info
	   set paramType [string tolower [lindex $paramSet 0]]
	   set paramName [lindex $paramSet 1]
	   set paramSet [lrange $paramSet 2 end]

#Extract dimension info.
	   set size ""
	   while {1} {
		if {[llength $paramSet] == 0} break
		set next [lindex $paramSet 0]
		if {[string index $next 0] != "<"} break
		if {![regexp {<([0-9]+)>} $paramSet all length]} break
		set paramSet [lrange $paramSet 1 end]
		lappend size $length
		}

#Size should be an int, even for scalar variables.
	   if {$size == ""} {set size 1}
	   lappend table($typeName) $paramName
	   set table($typeName,$paramName,type) $paramType
	   set table($typeName,$paramName,size) $size

#Extract constraint info
	   set constraint [lindex $paramSet 0]
	   set paramSet [lrange $paramSet 1 end]
	   set table($typeName,$paramName,constraint) $constraint

#Enum.  I could also add code to handle Dervish enum datatypes.
	   if {$constraint == "enum"} {
		set table($typeName,$paramName,enumlist) $paramSet
		}

#Add Dervish enum constraint
	   if {[lsearch $enum() $paramType] >= 0} {
		set table($typeName,$paramName,type) char
		set table($typeName,$paramName,size) 80
		set table($typeName,$paramName,constraint) enum
		set table($typeName,$paramName,enumlist) $enum($paramType)
		}

#Oneof.
	   if {$constraint == "oneof"} {
		set oneof [split [lindex $paramSet 0] :]
		set table($typeName,$paramName,oneoftable) [lindex $oneof 0]
		set table($typeName,$paramName,oneofattr) [lindex $oneof 1]
		}

#Default
	   if {$constraint == "default"} {
		set default [lindex $paramSet 0]
		set table($typeName,$paramName,default) $default
		}
#Other constraints don't need any extra processing
	   }
	set table($typeName,counter) 0
   } else {
	close $fid
	error "Bad typedef type: $typeType"
	}
   return
   }

#########################################################################
#Read and decode a parameter file.  Inputs are:
#	paramFile	Name of parameter file
#	hdrPtr		Scalar variable to contain header keyed list
#	tablePtr	Array variable to contain table defs and data
#	enumPtr		Array variable to contain enum defs (optional)

proc paramRead {paramFile hdrPtr tablePtr {enumPtr ""}} {

   set type NULL

   set structList ""
#
# setup the hdr variable
#

    upvar $hdrPtr hdr
    upvar $tablePtr table
    if {$enumPtr != ""} {upvar $enumPtr enum}
    foreach var "hdr table enum" {
	if {[info exists $var]} {unset $var}
	}

    set hdr ""
    set table() ""
    set enum() ""
#
# open the ascii file
#
    set fid [open $paramFile r]

#
# go through it line by line
#
   set lines ""
   set linebuf ""
   while {1} {
	set nextline [string trim [gets $fid]]
	if {[eof $fid]} break
#If we have a continuation at end of this line, stash what we have so far
#and fetch the next bit.  We pad end of line with a blank.
	if {[regexp {(.+)\\$} $nextline all nextline]} {
	   append linebuf "$nextline "
	   continue
	} else {
	   append linebuf $nextline
	   }
	set line $linebuf
	set linebuf ""

#Check for comment lines and blank lines
	if {[string index $line 0] == "#"} continue
	if {[string compare $line ""] == 0} continue

#Delete trailing comments
	regsub {#.*$} $line "" line

#If line is not complete, read more	
	append lines $line
	if {![info complete $lines]} {
#	   append lines " "
	   continue
	   }

#Reset multiple line buffer
	set line $lines
	set lines ""

#Convert multiple tabs and spaces to single blanks
	regsub -all "\t+" $line { } line
	regsub -all " +"  $line { } line

#Key is lower case for table records, typedef.
#Is this a typedef?
	set key [lindex $line 0]
	set keylc [string tolower $key]
	if {$keylc == "typedef"} {
	   paramReadTypedef $line $fid

#Is this line a record in a table?
	} elseif {[lsearch $table() $keylc] >= 0} {
	   set table($keylc,$table($keylc,counter)) [lrange $line 1 end]
	   incr table($keylc,counter)
	} else {

#Line is keyword-value pair for header.  Use real key name, not lower case version
	   set line [lrange $line 1 end]
	   keylset hdr $key $line
	   }
	}

   close $fid
#Anything left over:
   if {[string compare $linebuf ""] != 0} {
	error "Incomplete line found: $linebuf."
	}
   if {[string compare $lines ""] != 0} {
	error "Incomplete lines found: $lines."
	}
   paramParse table
   return
   }

#####################################################################
#Parse a table - make it easier to deal with.

proc paramParse {_intable} {
   upvar 1 $_intable table
   foreach type $table() {
	set names $table($type)
	loop j 0 [llength $names] {
	   set name [lindex $names $j]
	   loop i 0 $table($type,counter) {
		set val [lindex $table($type,$i) $j]
		set table($type,$i,$name) $val
		}
	   }
	}
   return
   }

