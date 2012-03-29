# Read and write FITS table files

################################################################################
# Routines for reading a FITS Table file

#Get size of an array.  Return 0 if it does not exist
proc arrsize {var} {
   upvar 1 $var local
   if {[catch {set size [array size local]}]} then {return 0} else {
                return $size}
   }

########################################################################
# Pad an input string with blanks to size characters.  Truncate if over.

proc pad {instring size} {
   set len [expr {$size - [string length $instring]}]
   for {set i 1} {$i <= $len} {incr i} {
        append instring " "
        }
   return [string range $instring 0 [expr {$size - 1}]]
   }

#######################################################################
proc pad80 {instring} {
   return [pad $instring 80]
   }

########################################################################
#Given a header, find the first line no. with a given keyword
proc HdrGetLineNo {inheader keyword} {
   upvar 1 $inheader header
   foreach i [lsort -integer [array names header]] {
	if {[HdrKeyGet $header($i)] == $keyword} then {return $i}
	}
   return 0
   }

########################################################################
#Insert pre-formatted line into header.  Insert at end - screw the
#concept of specifying the line number.
proc HdrInsLine {inheader line} {
   upvar 1 $inheader header
   if {[info exists header]} {
	set indices [array names header]
	set max [eval max $indices]
	set n [expr $max+1]
   } else {
	set n 1
	}
   set header($n) [pad80 $line]
   return
   }

########################################################################
proc hdrInsertLine {hdr n line} {
   set name $hdr
   global $name
#   if {![info exists $name]} {error "Can't find header $hdr"}
   HdrInsLine $hdr $line
   return
   }

########################################################################
#Extract the keyword field
proc HdrKeyGet {line} {
   return [lindex [string range $line 0 7] 0]
   }

########################################################################
#Extract the parameter type.  This is REEEALLL simple - N or A
proc HdrParTypeGet {line} {
   set line [string range $line 10 end]
#Get rid of leading and trailing spaces.
   set line [string trim $line " "]
   if {[string index $line 0] != {'} } then {
	return I
   } else {
	return A
	}
   }

########################################################################
#Extract the parameter field
proc HdrParGet {line} {
   set line [string range $line 10 end]
#Get rid of leading and trailing spaces.
   set line [string trim $line " "]
   if {[string index $line 0] != {'} } then {
	return [lindex $line 0]
	}
   set line [string range $line 1 end]
   set len [string length $line]
   set outline {}
   set flag 0
   for {set i 0} {$i < $len} {incr i} {
	set char [string index $line $i]
	if {$char == {"} } then {
	   append outline {\"}
	   continue
	   }
	if {$char == {'}} then {
	   if {$flag == 1} then {
		append outline {'}
		set flag 0
		continue
		}
	   set flag 1
	   } else {
	   if {$flag == 1} then {
		break}
	   append outline $char
	   continue
	   }
	}
   if {$flag == 1} then {return $outline} else  {
# Else error. For now, return the rest of the line as if the last quote were
# present
	return $outline}
   }

########################################################################
#Extract the comment field
proc HdrComGet {line} {
   global VERBOSE
   set line [string range $line 10 end]
#Fits recommends but does not mandate that character strings begin in col.11
#Hence the need to trim blanks
   set line [string trim $line " "]
   if {[string index $line 0] != {'} } then {
	lvarpop line
	set line [string trimleft $line / ]
	return $line
	}
   set line [string range $line 1 end]
   set len [string length $line]
   set flag 0
   for {set i 0} {$i < $len} {incr i} {
	set char [string index $line $i]
	if {$char == {'}} then {
	   if {$flag == 1} then {
		set flag 0
		continue
		}
	   set flag 1
	   } else {
	   if {$flag == 1} then {
		break}
	   continue
	   }
	}
   if {$flag == 1} then {
#Look for leading /
	set line [string range $line $i end]
	if {$line == ""} return
	set size [string length $line]
	for {set i 0} {$i < $size} {incr i} {
	   set char [string index $line $i]
	   if {$char == " "} continue
	   if {$char == "/"} then {
		set retstring [string range $line [expr {$i+1}] end]
		return [string trim $retstring " "]
	   } else {
		set retstring [string range $line $i end]
		return [string trim $retstring " "]
		}
	   }
	return ""
# Else error. For now, return a null string as if the last quote were present
   } else {
	if {[info exists VERBOSE]} then {echo  Error in decoding comment}
	return ""
	}
   }

########################################################################
#Read a header from a FITS file, given a file handle
proc HdrRead {fhandle inheader} {
   upvar 1 $inheader header
   if {[info exists header]} {unset header}
   set nline 0
   while {1} {
	for {set i 0} {$i < 36} {incr i} {
	   catch {unset line}
	   set line [read $fhandle 80]
	   if {[eof $fhandle]} {
		close $fhandle
		error "Unexpected eof in HdrRead"
		}
	   if {[string length [string trim $line]] == 0} continue
	   if {[HdrKeyGet $line] == "COMMENT"} then {continue}
	   if {[HdrKeyGet $line] == "HISTORY"} then {continue}
	   incr nline
	   set header($nline) $line
	   }
	if {[HdrGetLineNo header END] > 0} then  {return}
	}
   }

########################################################################
#Extract an ascii field from a header, given a keyword
proc HdrGetAsAscii {inheader keyword} {
   upvar 1 $inheader header
   foreach i [array names header] {
	if {[HdrKeyGet $header($i)] != $keyword} then {continue}
	return [HdrParGet $header($i)]
	}
   return ""
   }	

########################################################################
#Emulate dervish command.  hdr should be of format, eg., h0.hdr, i.e. reg.hdr
proc hdrGetAsAscii {hdr keyword} {
   set name $hdr
   global $name
   if {![info exists $name]} {error "Can't find header $hdr"}
   return [HdrGetAsAscii $hdr $keyword]
   }

########################################################################
proc HdrInsWithAscii {inheader keyword value {comment ""}} {
   upvar 1 $inheader header
   if {[info exists header]} {
	set indices [array names header]
	set max [eval max $indices]
	set n [expr $max+1]
   } else {
	set n 1
	}
   regsub -all ' $value '' value
   set value [format '%-8s' $value]
   if {[string length $value] > 60} {
	error "Value $value too big for FITS header field!"
	}
   set header($n) [format "%-8s= %-20s / %-47s" [string range $keyword 0 7] \
	$value [string range $comment 0 46]]
   set header($n) [pad80 $header($n)]
   return
   }

########################################################################
proc hdrInsWithAscii {hdr keyword value {comment ""}} {
   set name $hdr
   global $name
#   if {![info exists $name]} {error "Can't find header $hdr"}
   HdrInsWithAscii $name $keyword $value $comment
   return
   }

########################################################################
#Extract an integer field from a header, given a keyword
proc HdrGetAsInt {inheader keyword} {
   upvar 1 $inheader header
   foreach i [array names header] {
	if {[HdrKeyGet $header($i)] != $keyword} then {continue}
	return [HdrParGet $header($i)]
	}
   return ""
   }	

########################################################################
#Emulate dervish command.  hdr should be of format, eg., h0.hdr, i.e. reg.hdr
proc hdrGetAsInt {hdr keyword} {
   set name $hdr
   global $name
   if {![info exists $name]} {error "Can't find header $hdr"}
   return [HdrGetAsInt $hdr $keyword]
   }

########################################################################
proc HdrInsWithInt {inheader keyword value {comment ""}} {
   upvar 1 $inheader header
   if {[info exists header]} {
	set indices [array names header]
	set max [eval max $indices]
	set n [expr $max+1]
   } else {
	set n 1
	}
   if {[string length $value] > 60} {
	error "Value $value too big for FITS header field!"
	}
   if {[catch {format %d $value}]} {
	error "Value $value not an integer!"
	}
   set header($n) [format "%-8s= %20s / %-47s" [string range $keyword 0 7] \
	$value [string range $comment 0 46]]
   set header($n) [pad80 $header($n)]
   return
   }

########################################################################
proc hdrInsWithInt {hdr keyword value {comment ""}} {
   set name $hdr
   global $name
#   if {![info exists $name]} {error "Can't find header $hdr"}
   HdrInsWithInt $name $keyword $value $comment
   return
   }

########################################################################
#Extract a double field from a header, given a keyword
proc HdrGetAsDbl {inheader keyword} {
   upvar 1 $inheader header
   foreach i [array names header] {
	if {[HdrKeyGet $header($i)] != $keyword} then {continue}
	return [HdrParGet $header($i)]
	}
   return ""
   }	

########################################################################
#Emulate dervish command.  hdr should be of format, eg., h0.hdr, i.e. reg.hdr
proc hdrGetAsDbl {hdr keyword} {
   set name $hdr
   global $name
   if {![info exists $name]} {error "Can't find header $hdr"}
   return [HdrGetAsDbl $hdr $keyword]
   }

########################################################################
proc HdrInsWithDbl {inheader keyword value {comment ""}} {
   upvar 1 $inheader header
   if {[info exists header]} {
	set indices [array names header]
	set max [eval max $indices]
	set n [expr $max+1]
   } else {
	set n 1
	}
   if {[string length $value] > 60} {
	error "Value $value too big for FITS header field!"
	}
   if {[catch {format %f $value}]} {
	error "Value $value not a floating point number!"
	}
   set header($n) [format "%-8s= %20s / %-47s" [string range $keyword 0 7] \
	$value [string range $comment 0 46]]
   set header($n) [pad80 $header($n)]
   return
   }

########################################################################
proc hdrInsWithDbl {hdr keyword value {comment ""}} {
   set name $hdr
   global $name
#   if {![info exists $name]} {error "Can't find header $hdr"}
   HdrInsWithDbl $name $keyword $value $comment
   return
   }

########################################################################
#Extract a logical field from a header, given a keyword
proc HdrGetAsLogical {inheader keyword} {
   upvar 1 $inheader header
   foreach i [array names header] {
	if {[HdrKeyGet $header($i)] != $keyword} then {continue}
	return [HdrParGet $header($i)]
	}
   return ""
   }	

########################################################################
#Emulate dervish command.  hdr should be of format, eg., h0.hdr, i.e. reg.hdr
proc hdrGetAsLogical {hdr keyword} {
   set name $hdr
   global $name
   if {![info exists $name]} {error "Can't find header $hdr"}
   return [HdrGetAsLogical $hdr $keyword]
   }

########################################################################
proc HdrInsWithLogical {inheader keyword value {comment ""}} {
   upvar 1 $inheader header
   if {[info exists header]} {
	set indices [array names header]
	set max [eval max $indices]
	set n [expr $max+1]
   } else {
	set n 1
	}
   set value [string toupper $value]
   if {$value != "T" && $value != "F"} {
	error "Value $value is not a logical type!"
	}
   set header($n) [format "%-8s= %20s / %-47s" [string range $keyword 0 7] \
	[string range $value 0 19] [string range $comment 0 46]]
   set header($n) [pad80 $header($n)]
   return
   }

########################################################################
proc hdrInsWithLogical {hdr keyword value {comment ""}} {
   set name $hdr
   global $name
#   if {![info exists $name]} {error "Can't find header $hdr"}
   HdrInsWithLogical $name $keyword $value $comment
   return
   }

########################################################################
#Replace parameter value for a given header.  Keep original comment if
#no new one supplied.

proc HdrReplaceWithAscii {inheader keyword value {comment ""}} {
   upvar 1 $inheader header
   set index [HdrGetLineNo header $keyword]
   if {$index > 0} {
	set line $header($index)
	if {$comment == ""} {set comment [HdrComGet $line]}
	unset header($index)
	}
   HdrInsWithAscii header $keyword $value $comment
   return
   }

########################################################################
proc hdrReplaceWithAscii {hdr keyword value {comment ""}} {
   set name $hdr
   global $name
#   if {![info exists $name]} {error "Can't find header $hdr"}
   HdrReplaceWithAscii $name $keyword $value $comment
   return
   }

########################################################################
#Replace parameter value for a given header.  Keep original comment if
#no new one supplied.

proc HdrReplaceWithInt {inheader keyword value {comment ""}} {
   upvar 1 $inheader header
   set index [HdrGetLineNo header $keyword]
   if {$index > 0} {
	set line $header($index)
	if {$comment == ""} {set comment [HdrComGet $line]}
	unset header($index)
	}
   HdrInsWithInt header $keyword $value $comment
   return
   }

########################################################################
proc hdrReplaceWithInt {hdr keyword value {comment ""}} {
   set name $hdr
   global $name
#   if {![info exists $name]} {error "Can't find header $hdr"}
   HdrReplaceWithInt $name $keyword $value $comment
   return
   }

########################################################################
#Replace parameter value for a given header.  Keep original comment if
#no new one supplied.

proc HdrReplaceWithDbl {inheader keyword value {comment ""}} {
   upvar 1 $inheader header
   set index [HdrGetLineNo header $keyword]
   if {$index > 0} {
	set line $header($index)
	if {$comment == ""} {set comment [HdrComGet $line]}
	unset header($index)
	}
   HdrInsWithDbl header $keyword $value $comment
   return
   }

########################################################################
proc hdrReplaceWithDbl {hdr keyword value {comment ""}} {
   set name $hdr
   global $name
#   if {![info exists $name]} {error "Can't find header $hdr"}
   HdrReplaceWithDbl $name $keyword $value $comment
   return
   }

########################################################################
#Replace parameter value for a given header.  Keep original comment if
#no new one supplied.

proc HdrReplaceWithLogical {inheader keyword value {comment ""}} {
   upvar 1 $inheader header
   set index [HdrGetLineNo header $keyword]
   if {$index > 0} {
	set line $header($index)
	if {$comment == ""} {set comment [HdrComGet $line]}
	unset header($index)
	}
   HdrInsWithLogical header $keyword $value $comment
   return
   }

########################################################################
proc hdrReplaceWithLogcial {hdr keyword value {comment ""}} {
   set name $hdr
   global $name
#   if {![info exists $name]} {error "Can't find header $hdr"}
   HdrReplaceWithLogical $name $keyword $value $comment
   return
   }

########################################################################
#Add blank lines to a header to pad out to a 2880 byte record
proc HdrFill {inheader} {
   upvar 1 $inheader header
   set nline [arrsize header]
   while { [expr {($nline/36)*36 != $nline}]} {
        incr nline
        set header($nline) [pad80 ""]
        }
   }

########################################################################
proc hdrPrint {header} {
   global $header
   if {![info exists $header]} {
	error "No header $header"
	}
   foreach i [lsort -integer [array names $header]] {
	echo [string range [set ${header}($i)] 0 78]
	}
   return
   }

########################################################################
#### Read header and return a multitude of decoded header information. ###

proc fitsPduHdrDecode {inheader} {
   global VERBOSE
   upvar 1 $inheader header
   upvar 1 NROW NROW
   upvar 1 NCOL NCOL
   upvar 1 NBYTE NBYTE
   upvar 1 BITPIX BITPIX
   upvar 1 PCOUNT PCOUNT
   upvar 1 EXTEND EXTEND
   upvar 1 BSCALE BSCALE
   upvar 1 BZERO BZERO

#Init
   set PCOUNT 0
   set EXTEND 1
   if {[HdrGetAsLogical header SIMPLE] != "T"} then {
	if {[info exists VERBOSE]} then {
	   echo  "Not a SIMPLE=T file; proceeding anyway"
	   }
	}
   if {[HdrGetAsLogical header EXTEND] != "T"} then {
	set EXTEND 0
	}
   set BITPIX [HdrGetAsInt header BITPIX]
   set NBYTE [int [fabs [expr [HdrGetAsInt header BITPIX]/8]]]
   set NAXIS [HdrGetAsInt header NAXIS]
   if {  $NAXIS > 2} then {
	error "Too many axes"
	}
   if {$NAXIS >= 0} then {
	set NROW 0
	set NCOL 0
	}
   if {$NAXIS >= 1} then {
	set NROW 1
	set NCOL [HdrGetAsInt header NAXIS1]
	}
   if {$NAXIS >= 2} then {
	set NCOL [HdrGetAsInt header NAXIS1]
	set NROW [HdrGetAsInt header NAXIS2]
	}
   set BSCALE [HdrGetAsDbl header BSCALE]
   set BZERO [HdrGetAsDbl header BZERO]
   if {$BSCALE == ""} {set BSCALE 1}
   if {$BZERO == ""} {set BZERO 0}

#Delete all reserved lines
   set keywords {SIMPLE XTENSION BITPIX NAXIS NAXIS1 NAXIS2 PCOUNT \
	GCOUNT BZERO BSCALE END}
   set zaplist ""
   foreach keyword $keywords {
	lappend zaplist [HdrGetLineNo header $keyword]
	}
   foreach zap $zaplist {
	if {![info exists header($zap)]} continue
	unset header($zap)
	}
   return
   }

########################################################################
#Print out non-reserved header keywords
proc fitsTableHdrPrint {inheader} {

#The following variables are passed back to the calling routine by reference.
#Too many to pass as dummy arguments, so the calling program must be prepared
#to accept their absolute names
   global VERBOSE
   upvar 1 $inheader header
#Delete all reserved lines, then print the reset
   set nfield [HdrGetAsInt header TFIELDS]
   set keywords {XTENSION BITPIX NAXIS NAXIS1 NAXIS2 PCOUNT GCOUNT TFIELDS}
   set zaplist ""
   foreach keyword $keywords {
	lappend zaplist [HdrGetLineNo header $keyword]
	}
   loop k 1 [expr $nfield+1] {
	lappend zaplist [HdrGetLineNo header TFORM$k]
	lappend zaplist [HdrGetLineNo header TTYPE$k]
	lappend zaplist [HdrGetLineNo header TDIM$k]
	}
   foreach zap $zaplist {
	catch {unset header($zap)}
	}
   set nums [array names header]
   set nums [lsort -integer $nums]
   foreach num $nums {
	if {[string trim $header($num)] == ""} continue
	echo "   [string trim $header($num)]"
	}
   }

################################################################
#Read a FITS file - create region on the fly.

proc regReadFromFits {file} {
   if {![file exists $file]} {error "File $file does not exist!"}
   set fid [open $file]
   HdrRead $fid header
   fitsPduHdrDecode header
   set here [tell $fid]

#Need to pass the current file position to regReadFromFile because
#Tcl_GetOpenFile, at least on Linux, seems to mess with the actual
#position.

   set reg [regReadFromFile $fid $NROW $NCOL $BITPIX $BSCALE $BZERO $here]
   close $fid

#Final name for header array
#Stroke of genius!  Array names will be of form h1.hdr where h1 is region
#handle!
   set list [array get header]
   set name $reg.hdr
   global $name
   array set $name $list

#Set name
   set fname [string range [file tail $file] 0 78]
   handleSet $reg.name $fname
   return $reg
   }

################################################################
#Delete a region

proc regNew {nrow ncol} {
   set reg [regNewToReg $nrow $ncol]
   set name $reg.hdr
   global $name
   if {[info exists $name]} {unset $name}
   handleSet $reg.name $reg
   return $reg
   }

################################################################
#Delete a region

proc regDel {reg} {
   regDelFromReg $reg
   set name $reg.hdr
   global $name
   if {[info exists $name]} {unset $name}
   return
   }

################################################################
#Recreation of regReadFromFits, but create region on the fly always.

proc regWriteAsFits {reg file} {
   regWriteToFits $reg $file
   return
   }

#########################################################################
proc regWriteToFits {reg file} {
   set name $reg.hdr
   global $name

#Recreate important keywords
   HdrInsWithLogical outhead SIMPLE T
   HdrInsWithInt outhead BITPIX -32 "Bits per pixel"
   HdrInsWithInt outhead NAXIS 2 "Number of Axes"

#OK, let's get this right.  In a FITS file, the first axis increments first.
#This is column-direction.
   HdrInsWithInt outhead NAXIS1 [exprGet $reg.ncol] "Number of cols"
   HdrInsWithInt outhead NAXIS2 [exprGet $reg.nrow] "Number of rows"
   if {[info exists $name]} {
	if {[llength [array names $name]] > 0} {
	   foreach i [lsort -integer [array names $name]] {
		HdrInsLine outhead [set ${name}($i)]
		}
	   }
	}
   HdrInsLine outhead END
   HdrFill outhead
   set fid [open $file w]
   foreach i [lsort -integer [array names outhead]] {
	puts $fid $outhead($i) nonewline
	}
   flush $fid
   set location [tell $fid]
   regWrite $reg $fid
   close $fid
   return
   }

#######################################################################
#Emulate Dervish Command
proc subRegNew {reg nrow ncol row0 col0} {
   set name1 $reg.hdr
   global $name1
   set sreg [subRegFromReg $reg [expr round($row0)] [expr round($col0)] \
	   [expr round($nrow)] [expr round($ncol)]]
   set name2 $sreg.hdr
   global $name2
   set line [array get $name1]
   array set $name2 $line
   return $sreg
   }
#######################################################################
proc regList {} {
   echo [handleListFromType REGION]
   return
   }

#######################################################################
#Compress a region and copy the header along at the same time
proc compress {reg1 factor} {
   set name1 $reg1.hdr
   global $name1
   set reg2 [regCompress $reg1 $factor]
   set name2 $reg2.hdr
   global $name2
   array set $name2 [array get $name1]

#Need to change lots of header info
#   echo Need to change lots of header info!
   set scale [HdrGetAsDbl $name2 SCALE]
   if {$scale != ""} {
	set scale [expr $scale*$factor]
	HdrReplaceWithDbl $name2 SCALE $scale
	}
   set pixmm [HdrGetAsDbl $name2 PIXMM]
   if {$pixmm != ""} {
	set pixmm [expr $pixmm*$factor]
	HdrReplaceWithDbl $name2 PIXMM $pixmm
	}
   set flux [HdrGetAsDbl $name2 FLUX]
   if {$flux != ""} {
	set flux [expr $flux/($factor*$factor)]
	HdrReplaceWithDbl $name2 FLUX $flux
	}
   return $reg2
   }

#######################################################################

proc status {} {
   set list [handleListFromType REGION]
   echo [format "%5s%1s %-20s %8s %8s" ID " " Name Nrow Ncol]
   foreach hndl $list {
	set name [lindex [exprGet $hndl.name] 0]
	set nrow [exprGet $hndl.nrow]
	set ncol [exprGet $hndl.ncol]
	echo [format "%5s  %-20s %8d %8d" $hndl \
	   $name $nrow $ncol]
      }
   puts stdout "Memory in use: [memBytesInUse] bytes"
   puts stdout "Number of handles: [llength [handleList]]"
   return
   }

####################################################################
#Fetch a single pixel value

proc regPixGet {reg row col} {
   set nrow [exprGet $reg.nrow]
   set ncol [exprGet $reg.ncol]
   if {$row < 0 || $row >= $nrow} {error "row $row outside boundaries"}
   if {$col < 0 || $col >= $ncol} {error "col $col outside boundaries"}
   set n [expr $row*$ncol + $col]
   return [exprGet $reg.pixels<$n>]
   }

####################################################################
#set a single pixel value

proc regPixSet {reg row col val} {
   set nrow [exprGet $reg.nrow]
   set ncol [exprGet $reg.ncol]
   if {$row < 0 || $row >= $nrow} {error "row $row outside boundaries"}
   if {$col < 0 || $col >= $ncol} {error "col $col outside boundaries"}
   set n [expr $row*$ncol + $col]
   handleSet $reg.pixels<$n> $val
   return
   }

#########################################################################

#We don't have regCopy - so let's embed one region in another.
#in is input region to embed.
#out is output region in which we do the embedding
#irow0 and icol0 are the starting corner coords in the "out" region.

proc regEmbed {in out irow0 icol0} {
   set nrin [exprGet $in.nrow]
   set ncin [exprGet $in.ncol]
   set nrout [exprGet $out.nrow]
   set ncout [exprGet $out.ncol]
   loop i 0 $nrin {
	loop j 0 $ncin {
	   set irow [expr $i + $irow0]
	   if {$irow < 0 || $irow >= $nrout} continue
	   set icol [expr $j + $icol0]
	   if {$icol < 0 || $icol >= $ncout} continue
	   regPixSet $out $irow $icol [regPixGet $in $i $j]
	   }
	}
   return
   }
