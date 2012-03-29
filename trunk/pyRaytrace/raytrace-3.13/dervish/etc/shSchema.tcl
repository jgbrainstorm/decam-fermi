# schemaPrint <handle>	Print the schema for a handle                         #
# schemaPrintFromType <type>                                                  #
#                       Print the schema for a type                           #
# schemaMemberPrint <type> <member>                                           #
#                       Print the schema for the member of a type             #
###############################################################################
#
# print a type's schema
#
proc schemaPrint {handle} {
	schemaPrintFromType [handleType $handle]
}

###############################################################################
proc schemaPrintFromType {type} {
	foreach f [schemaGetFromType $type] { echo $f }
}

###############################################################################
proc schemaMemberGet {type elem} {
	set members [schemaGetFromType $type]
	keylget members $elem
}

#Procedures defined here:
#	schemaRecurseGetFromType type
#	recursTrans inresult arg type nstar dim offset  (private routine)
#	schemaRecursePrintFromType type
#	typeSizeOf
#	typePad
#	typedef
#	makeio
#	makeioSub
#       schemaGetFull
###############################################################################

proc schemaGetFull { handle } {
    return [schemaGetFullFromType [handleType $handle]]
}

set fits2SchemaArgs {
    {fits2Schema "Read in a FITS table file and create the appropriate objects.\n" }
    {<file> STRING "" file "FITS file to be read"}
    {<type> STRING "" type "name of schema to create"}
    {[hdr]  STRING NULL hdr "optional header, set if not NULL"}
    {-hdu   STRING 1 hdu "hdu to read in file"}
    {-sc   STRING "" sc "use this schemaTrans instead of auto-generating"}    
    {-existingType   CONSTANT "" existingType "use an existing type"}
}
ftclHelpDefine shFitsIo fits2Schema [shTclGetHelpInfo fits2Schema $fits2SchemaArgs]

#############################################################################

#schemaGet routines
##############################################################################
#
proc schemaRecurseGetFromType {type} {
   set result ""
   recursTrans result "" $type 0 NULL ""
   return $result
   }

#######################################################################

proc recursTrans {inresult arg type nstar dim offset} {
upvar 1 $inresult result
   set test [schemaKindGetFromType $type]
   if {$test == "PRIM" || $test == "ENUM" || $nstar != 0} {
# type is primitive.
	lappend result [list $arg $type $nstar $dim $offset]
	return
	}
# type is complex.  Is it an array?  If so, expand each array element.
   if {$dim != "NULL"} then {
	set ndim [llength $dim]
	set dimtot 1
	for {set i [expr $ndim-1]} {$i >= 0} {set i [expr $i-1]} {
	   set dim1($i) $dimtot
	   set dim2($i) [lindex $dim $i]
	   set dimtot [expr $dimtot*$dim2($i)]
	   }
	loop j 0 $dimtot {
	   set subarg $arg
	   for {set i 0} {$i < $ndim} {incr i} {
		set subdim [fmod [floor [expr $j/$dim1($i)]] $dim2($i)]
		append subarg <$subdim>
		}
	   recursTrans result $subarg $type $nstar NULL $offset
	   }
	return
	}	   
#Complex inline type.	
   set allschema [schemaGetFullFromType $type]
   lvarpop allschema
   foreach f $allschema {
	set subarg [lindex $f 0]
	set subtype [lindex $f 1]
	set substar [lindex $f 2]
	set subdim [lindex $f 3]
	set suboff [lindex $f 4]
	recursTrans result "$arg $subarg" $subtype $substar $subdim $suboff
	}
   }

#######################################################################

proc schemaRecursePrintFromType {type} {
   set biglist [schemaRecurseGetFromType $type]
   set nelem [llength $biglist]
   for {set i 0} {$i < $nelem} {incr i} {
	set list [lindex $biglist $i]
	set arglist [lindex $list 0]
	set type [lindex $list 1]
	set nstar [lindex $list 2]
	set dim [lindex $list 3]
	set offset [lindex $list 4]
	set narg [llength $arglist]
	set arg [lindex $arglist 0]
	set j 1
	while {$j < $narg} {
	   set arg $arg.[lindex $arglist $j]
	   incr j
	   }
	set j 0
	while {$j < $nstar} {
	   append type *
	   incr j
	   }
	if {$dim != "NULL"} then {
	   for {set j 0} {$j < [llength $dim]} {incr j} {
		set type $type\[[lindex $dim $j]\]
		}
	   }
	echo [format "%-30s%s" $arg $type]
	}
   }


##########################################################################
#More easily used interface to schemaDefine facility.  Permit a more C-like
#syntax for defining dervish TYPEs at run-time
#
#example:
#   typedef struct {
#	int i;
#	float a<12>;
#	REGION *reg;
#	MASK mask;
#	} MYTYPE;
#
#Restrictions: Embedded structs like mask above must be predefined.
#Arrays can be specified either with [] or <> notation
#
# Enums should work!
#
#For a given data type, find the maximum size of each embedded primitive.
#C compiler will pad to the largest embedded primitive type!

############################################################################
proc typeSizeOf {type} {
   set schema [schemaGetFullFromType $type]
   set schemaElem [lindex $schema 0]
   set size [lindex $schemaElem 2]
   return $size
   }

############################################################################
proc typepad {type} {
   set schema [schemaRecurseGetFromType $type]
   set max 0
   loop i 0 [llength $schema] {
	set schemaElem [lindex $schema $i]
	set schemaType [lindex $schemaElem 1]
	set nstar [lindex $schemaElem 2]
	if {$nstar > 0} then {set schemaType PTR}
	set max [max $max [typeSizeOf $schemaType]]
	}
   return $max
   }
   
############################################################################
proc typedef {args} {
   set nsch 1
#There are different valid formats for typedef's.
   set len [llength $args]
#2 arguments: example - typedef SECTION RANGE.  This cannot be represented
#in dervish.  The TCL proc makeio parses legal C structures and sends them
#here with possible only 2 args.  Thus we need to recover gracefully.
   if {$len == 2} {
	echo "Bypassing 2-argument input to typedef"
	return
	}
   if {$len != 3 && $len != 4} then {error "Wrong number of arguments"}
   if {$len == 3} then {
	set struct [lindex $args 0]
	set structdef [lindex $args 1]
	set type [lindex $args 2]
	set tag $type
	}
#If 4 args, we have a tag (not needed)
   if {$len == 4} then {
	set struct [lindex $args 0]
	set tag [lindex $args 1]
	set structdef [lindex $args 2]
	set type [lindex $args 3]
	}
   if {$struct != "struct" && $struct != "enum"} then {
	echo "Can only create STRUCTs or ENUMs at present"
	return
	}
# <cr>'s are likely part of the structdef argument.  Purge here
   regsub -all \n $structdef ";" structdef
   regsub -all \t $structdef " " structdef
   regsub -all \;\; $structdef \; structdef
   set pad 0
#Find and zonk all comments
   while {1} {
	if {![regexp -indices (/\[*\]) $structdef c1]} then break
	set first [string range $structdef 0 [expr [lindex $c1 0]-1]]
	set last [string range $structdef [expr [lindex $c1 1]+1] end]
	if {![regexp -indices (\[*\]/) $last c2]} then {error "unmatched /* "}
	set last [string range $last [expr [lindex $c2 1]+1] end]
	set structdef ${first}${last}
	}
   set offset 0
   set nelem 0
#Do enums here
   if {$struct == "enum"} then {
	set value 0
#Split the definition based on commas
	set split [split $structdef ,]
	loop i 0 [llength $split] {
	   set structElem [lindex $split $i]
	   regsub = $structElem " = " structElem
	   if {[llength $structElem] == 0} then continue
	   set name [lindex $structElem 0] 
#Do we have an = sign?
	   if {[lindex $structElem 1] == "="} then {
		set value [lindex $structElem 2]
#Allow for variable substitution
		catch {global $value}
		if {[info exists $value]} then {
		   set value [set $value]}
		}
	   if {[catch {format %d $value}]} then {
		echo Non-numerical value $value used for enum
		}
	   incr nelem
	   set schemaElem($nsch)($nelem) "$name $type $value NULL 0"
	   set $name $value
	   incr value
	   }
	set size [typeSizeOf int]
	set schema($nsch) "$type ENUM $size $nelem"
	schemaDefine schema schemaElem
	return
	}
#Split the definition based on semicolons
   set split [split $structdef \;]
   loop i 0 [llength $split] {
	set structElem [lindex $split $i]
#Save for error message later, want original, not one with substitutions. (efb)
        set localStructElem $structElem
        set localPad 0
	if {[llength $structElem] == 0} then continue
#A complex line might be   unsigned short int ** abc<4>;
#The ** pattern need not have spaces on either side.
#Normalize here
#Delete trailing spaces around a ** pattern
	regsub -all "(\[*]+)\[ ]*" $structElem \\1 structElem
#Insert leading spaces
	regsub -all "(\[^*])\[*]" $structElem "\\1 *" structElem
#Likewise for the dimension information
	regsub -all \[\[\] $structElem < structElem
	regsub -all \[\]\] $structElem > structElem
#Delete all leading space before the open <
	regsub -all "\[ ]*<\[ ]*" $structElem < structElem
	regsub -all "\[ ]*>\[ ]*" $structElem > structElem
#Replace short int with short, long int with long.  Curse these variations..
	regsub -all "short\[ ]+int" $structElem short structElem
	regsub -all "long\[ ]+int" $structElem short structElem
	set elemType [lvarpop structElem]
	if {$elemType == "const"} then {set elemType [lvarpop structElem]}
	if {$elemType == "unsigned"} then {
		set elemType [list $elemType [lvarpop structElem]]
		}
	if {$elemType == "signed"} then {
		set elemType [list $elemType [lvarpop structElem]]
		}
	if {$elemType == "struct"} then {
	   set elemType [lvarpop structElem]
	   if {$elemType == "$tag"} then {set elemType $type}
	   }
	set vars $structElem
#Might be multiple variables specified
	regsub -all , $vars " " vars
#Parse indirections
	loop i 0 [llength $vars] {
	   set var [lindex $vars $i]
	   set nstar 0
	   while {1} {
		if {[string index $var 0] != "*"} then break
		incr nstar
		set var [string range $var 1 end]
		}
#Parse array dimensions
	   set ndim 0
	   set dim ""
	   set dimtot 1
	   regsub -all ">" $var "" var
	   set var [split $var <]
	   set name [lindex $var 0]
	   loop j 1 [llength $var] {
		incr ndim
		set idim [lindex $var $j]
#Allow for variable substitutions
#Look up 1 level for a variable definition
#		upvar 1 $idim cvar
#Nah, make it global.  Too much nesting gets in the way otherwise
		catch {global $idim}
		if {[info exists $idim]} then {
		   set idim [set $idim]}
		set dimtot [expr $dimtot*$idim]
		append dim $idim " "
		}
	   set dim [string trimright $dim " "]
	   if {$ndim == 0} then {set dim NULL}
	   incr nelem
#Pad according to the current data type.
	   if {$nstar == 0} then {
		set thispad [typepad $elemType]
	   } else {
		set thispad [typepad PTR]
		}
	   set offset [int [expr ([floor [expr ($offset-1)/$thispad]]+1)*$thispad]]
	   set schemaElem($nsch)($nelem) [list $name $elemType $nstar $dim $offset]
#Get size of 1 element
	   if {$nstar == 0} then {
		set size [typeSizeOf $elemType]
	   } else {
		set size [typeSizeOf PTR]
		}
	   set totsize [expr $size*$dimtot]
	   set offset [expr $offset+$totsize]
	   if {$nstar == 0} then {
# Save pad for this variable in order to detect errors in variable
# syntax.
	        set localPad [typepad $elemType]
	   } else {
	        set localPad [typepad PTR]
		}
	   set pad [max $pad $localPad]
	   }
# if we did not get a value for localPad, then there was an error. this is
# probably an error in the expression syntax. (efb)
           if {$localPad == 0} then {
	       echo "Invalid expression - $localStructElem"
	       return
	   }
	}
   set offset [int [expr ([floor [expr ($offset-1)/$pad]]+1)*$pad]]
   set schema($nsch) "$type STRUCT $offset $nelem"
   # schemaDefine may return an error.  Since this is our last 
   # statement before an information-less return we can just let
   # tcl handle this error.
   schemaDefine schema schemaElem
   return
   }

##########################################################################
#
#Top level routines operates on a group of .h files and load schema in one go.
#This is how the dervish makeio facility operates.

proc makeio {args} {
#Check 1st argument.  If it is -pragma, remove from the argument list and
#pass it as a flag to makeioSub
   if {[lindex $args 0] == "-pragma"} then {
	lvarpop args
	set flag 1
   } else {
	set flag 0
	}
   foreach arg $args {
	set files [glob $arg]
	foreach file $files {
	   makeioSub $flag $file
	   }
	}
   return
   }

#####################################################################
#Read a .h file and parse off the typedef's.
#If flag is 1, we check for pragmaSCHEMA and relatives.
#
proc makeioSub {flag file} {
   set fid [open $file]
   while {1} {
	set line [gets $fid]
	if {[eof $fid]} then break
#C allows braces to not have surrounding whitespace - an anathema to TCL
#Fix that here.  Also, convert tabs to spaces
	regsub -all \[\{] $line " \{ " line
	regsub -all \[\}] $line " \} " line
	regsub -all \t $line " " line
#Should also take care of character strings here
#Search for #define statements.
	if {[regsub #define $line "" line]} then {
	   set var [lindex $line 0]
	   set val [lindex $line 1]
#Let us make the variable name global; make it available to other files
	   global $var
	   set $var $val
	   continue
	   }
	append list $line " "
	}
   close $fid
#Find and zonk all comments.  I might have a "typedef" embedded in one!
#Convert all /* pragma SCHEMA */ statements to #schema#, which is not valid C
#Convert all /* pragma CONSTRUCTOR */ statements to #schema#, which is not valid
#No constructors, however, are allowed.
   regsub -all \n $list " " list
   regsub -all "pragma\[ ]*NOCONSTRUCTOR" $list \
	" " list
   regsub -all "/\[*]\[ ]*pragma\[ ]*SCHEMA\[ ]*\[*]/" $list \
	"#schema#" list
   regsub -all "/\[*]\[ ]*pragma\[ ]*CONSTRUCTOR\[ ]*\[*]/" $list \
	"#schema#" list
   regsub -all "/\[*]\[ ]*pragma\[ ]*USER\[ ]*\[*]/" $list \
	"#schema#" list
   regsub -all "/\[*]\[ ]*pragma\[ ]*LOCK\[ ]*\[*]/" $list \
	"#schema#" list
   while {1} {
	if {![regexp -indices (/\[*\]) $list c1]} then break
	set first [string range $list 0 [expr [lindex $c1 0]-1]]
	set last [string range $list [expr [lindex $c1 1]+1] end]
	if {![regexp -indices (\[*\]/) $last c2]} then {error "unmatched /* "}
	set last [string range $last [expr [lindex $c2 1]+1] end]
	set list "${first} ${last}"
	}
#Find and zonk all character strings.
   while {1} {
#First, zonk any embedded quotes
	regsub -all \\\" $list " " list
	if {![regexp -indices \" $list c1]} then break
	set first [string range $list 0 [expr [lindex $c1 0]-1]]
	set last [string range $list [expr [lindex $c1 1]+1] end]
	if {![regexp -indices \" $last c2]} then {error "unmatched \" "}
	set last [string range $last [expr [lindex $c2 1]+1] end]
	set list "${first} ${last}"
	}
#search for a "typedef" string
   while {1} {
	if {![regexp -indices "typedef " $list c1]} then break
	set list [string range $list [lindex $c1 0] end]
#We could use lvarpop, but it is slow
	set command [lindex $list 0]
	set list [lrange $list 1 end]
#Because of possible variations on specification of typedefs, I will
#just look for the next list item that ends with a semicolon.
#We could use lvarpop, but it is slow.
	while {1} {
	   set item [lindex $list 0]
	   if {$item == ""} then continue
	   set list [lrange $list 1 end]
	   if {[regexp -indices (\;\$) $item c1]} then {
		set item [string trimright $item \;]
		if {$item != ""} then {lappend command $item}
		break
		}
	   lappend command $item
	   }
#If flag is 1, we check for a #schema# following the typedef
	if {$flag == 0 || [lindex $list 0] == "#schema#"} then {
echo Defining [lindex $command [expr [llength $command]-1]]
	   eval $command
	   }
	}
   return
   }
##########################################################################
ftclHelpDefine shSchema schemaRecurseGetFromType "schemaRecurseGetFromType type"
ftclHelpDefine shSchema schemaRecursePrintFromType \
   "schemaRecursePrintFromType type"
ftclHelpDefine shSchema typedef "\{struct | enum\} \[tag\] \{definition} type;"
ftclHelpDefine shSchema makeio "filespec \[filespec ...\]"
ftclHelpDefine shSchema schemaGetFull "Usage: schemaGetFull <handle expr>\n\t\t\
Return the complete definition of a schema in the form\n\t\t\
used by schemaDefine."
return
