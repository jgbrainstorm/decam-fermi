#Helper procedures to prepare glass info catalogs.
#Because different vendors use different formats for info in catalogs,
#everything must be hand crafted.

#Top commands:	glassTransPrep
#		glassDispPrep
#		glassDenPrep

global here
set here [file dirname [info script]]

########################################################################
#Read Zemax file.
#outid is the output file handle, _innames is a temporary global variable
#to store the ensemble of all glass names.

proc zemaxTransRead {file outid _innames} {
   upvar 1 $_innames names
   global glassName
   global here
   set format "%-12s = %s"
   set format2 "t %-10s = %s"
   set maker [string tolower [file root [file tail $file]]]
   regsub -all {[0-9]+} $maker "" maker
   regsub -all {_} $maker "" maker
   set fid [open $file]
   set nm ""
   set gwaves ""
   set gtrans ""
   while {1} {
	set line [gets $fid]
	if {[eof $fid]} break
	set code [lindex $line 0]
	if {$code == "NM"} {
	   if {[llength $gwaves] > 0 && [lsearch $names $glass] < 0} {
		puts $outid ""
		puts $outid [format $format name $glass]
		puts $outid [format $format maker schott]
		for {set i 0} {$i < [llength $gwaves]} {incr i} {
		   set wave [lindex $gwaves $i]
		   set trans [lindex $gtrans $i]
		   puts $outid [format $format2 $wave $trans]
		   }
		lappend names $glass
		}
	   set glass [lindex $line 1]
#For Hoya, insert "-" between letters and numbers
	   if {$maker == "hoya"} {
		regsub -all {([A-Z]+)([0-9]+)} $glass {\1-\2} glass
		}

#Convert to case sensitive name
	   if {[info exists glassName($glass)]} {
		set glass $glassName($glass)
	   } else {
		echo No name translation for $glass

#One fixup that I know about.
		regsub LASF $glass LaSF glass
		}
	   set gwaves ""
	   set gtrans ""
	   }
	if {$code == "IT"} {
	   set wave [lindex $line 1]
	   set trans [lindex $line 2]
	   set thick [lindex $line 3]

#Rescale everything.
	   set wave [expr $wave*1.e4]
	   set trans [format %.3f [expr pow($trans,$thick/25.)]]
	   lappend gwaves $wave
	   lappend gtrans $trans
	   lappend gthick $thick
	   }
	}
   close $fid


#Clean up last filter.
   if {[llength $gwaves] > 0} {
	puts $outid ""
	puts $outid [format $format name $glass]
	puts $outid [format $format maker schott]
	for {set i 0} {$i < [llength $gwaves]} {incr i} {
	   set wave [lindex $gwaves $i]
	   set trans [lindex $gtrans $i]
	   puts $outid [format $format2 $wave $trans]
	   }
	lappend names $glass
	}
   return
   }

#######################################################################
#Read Transmission data
proc glassTransPrep {} {
   global here glassName

#Get translation list
   echo fetching translation list
   opusName
   set outfile $here/TRANSMIT.TXT
   set outid [open $outfile w]
   set format "%-12s = %s"
   set format2 "t %-10s = %s"
   set names ""

#Original Schott - read from ZMAX file.
   set file $here/ZEMAX/SCHOTT.AGF
   zemaxTransRead $file $outid names

#Fake SIO2 - read from ZMAX file.
   set file $here/ZEMAX/SIO2.AGF
   zemaxTransRead $file $outid names

#New Schott.  I created a stripped-down .csv from the excel spreadsheet.
   set file $here/../doc/schott.trans.csv
   set fid [open $file]
   set line [split [gets $fid] ,]
   set n [llength $line]
   set waves dummy
   loop i 1 $n {
	set title [lindex [lindex $line $i] 0]
	regsub TAUI25/ $title "" title
	lappend waves $title
	}
   while {1} {
	set line [split [gets $fid] ,]
	if {[eof $fid]} break
	if {[llength $line] == 0} continue
	set glass [lindex [lindex $line 0] 0]
	if {[lsearch $names $glass] >= 0} continue
	lappend names $glass
	set gwaves ""
	set gtrans ""
	loop i 1 $n {
	   set wave [lindex $waves $i]
	   set trans [lindex $line $i]
	   if {$trans != ""} {
		lappend gwaves $wave
		lappend gtrans $trans
		}
	   }
	if {[llength $gwaves] > 0} {
	   puts $outid ""
	   puts $outid [format $format name $glass]
	   puts $outid [format $format maker schott]
	   for {set i [expr [llength $gwaves]-1]} {$i >= 0} {incr i -1} {
		set wave [lindex $gwaves $i]
		set trans [lindex $gtrans $i]
		set wave [expr $wave*10]
		puts $outid [format $format2 $wave $trans]
		}
	   }
	}
   close $fid

#New Schott - read from ZMAX file and catch any misses.
   set file $here/ZEMAX/SCHOTT_2000.AGF
   zemaxTransRead $file $outid names

#Ohara.  CSV file comes from stripped-down spreadsheet
   set file $here/../doc/ohara.trans.csv
   set fid [open $file]
   set line [split [gets $fid] ,]
   set n [llength $line]
   set waves "dummy dummy"
   loop i 2 $n {
	set title [lindex [lindex $line $i] 0]
	lappend waves $title
	}

#I read only the new glass name, but I should also process the old glass
#name, since this is what older designs will use.
   global glassName
   while {1} {
	set line [split [gets $fid] ,]
	if {[eof $fid]} break
	if {[llength $line] == 0} continue
	set glass [lindex [lindex $line 0] 0]
	if {[lsearch $names $glass] >= 0} continue
	regsub -all " " $glass "-" glass
	if {![info exists glassName($glass)]} {
	   echo No OPUS name for O'Hara $glass
	   }
	set gwaves ""
	set gtrans ""
	loop i 2 $n {
	   set wave [lindex $waves $i]
	   if {$wave > 1100} continue
	   set trans [lindex $line $i]
	   if {$trans != ""} {
		lappend gwaves $wave
		lappend gtrans $trans
		}
	   }
	if {[llength $gwaves] > 0} {
	   puts $outid ""
	   puts $outid [format $format name $glass]
	   puts $outid [format $format maker ohara]
	   for {set i 0} {$i < [llength $gwaves]} {incr i} {
		set wave [lindex $gwaves $i]
		set trans [lindex $gtrans $i]

#Tansmission is for 10 mm glass; change to 25 mm
		set trans [format %.3f [expr pow($trans,2.5)]]
		set wave [expr $wave*10]
		puts $outid [format $format2 $wave $trans]
		}
	   }
	}
   close $fid

#Newer O'Hara
   set file $here/ZEMAX/OHARA_2002.agf
   zemaxTransRead $file $outid names

#Hoya
   set file $here/ZEMAX/HOYA.AGF
   zemaxTransRead $file $outid names

#Sumita
   set file $here/ZEMAX/Sumita.agf
   zemaxTransRead $file $outid names

   close $outid
   return
   }

######################################################################
#Read OPUSGLASS.TXT and compile a list of names to be used to convert
#Zemax names to case-correct names.

proc opusName {} {
   global glassName
   if {[info exists glassName]} {unset glassName}
   global env
   set file $env(CRAY_DIR)/etc/OPUSGLASS.TXT
   set fid [open $file]
   set name ""
   while {1} {
	set line [string trim [gets $fid]]
	if {[eof $fid]} break
	if {[string length $line] == 0} continue
	set code [lindex $line 0]
	if {$code == "name"} {
	   set name [lindex $line 2]
	   continue
	   }
	if {$code == "maker"} {
	   set maker [lindex $line 2]
	   if {$name != ""} {
		set upper [string toupper $name]
		set glassName($upper) $name
		set glassName($name,maker) $maker
##		echo $upper: $name from $maker
		}
	   }
	}
   close $fid
   return
   }

##############################################################
########################################################################
#Read Zemax file.
#outid is the output file handle, _innames is a temporary global variable
#to store the ensemble of all glass names.

#This proc handles dispersion constants.
proc zemaxDispRead {file outid _innames} {
   upvar 1 $_innames names
   global glassName
   global here
   set format "%-8s %-8s %s %s"
   set maker [string tolower [file root [file tail $file]]]
   regsub -all {[0-9]+} $maker "" maker
   regsub -all {_} $maker "" maker
   set fid [open $file]
   set nm ""
   set type ""
   set coeffs ""
   while {1} {
	set line [gets $fid]
	if {[eof $fid]} break
	set code [lindex $line 0]
	if {$code == "NM"} {

#Output previous glass if we supplied name and coeffs.
	   if {[llength $coeffs] > 0 && [lsearch $names $glass] < 0} {
		puts $outid [format $format $glass $maker $type $coeffs]
		lappend names $glass
		}
	   set glass [lindex $line 1]

#For Hoya, insert "-" between letters and numbers
	   if {$maker == "hoya"} {
		regsub -all {([A-Z]+)([0-9]+)} $glass {\1-\2} glass
		}
	   set type [lindex $line 2]
	   if {$type == 1} {set type c} else {set type s}

#Convert to case sensitive name
	   if {[info exists glassName($glass)]} {
		set glass $glassName($glass)
	   } else {
		echo No name translation for $glass

#One fixup that I know about.
		regsub LASF $glass LaSF glass
		}
	   set coeffs ""
	   }
	if {$code == "CD"} {
	   loop i 1 7 {
		set coeff($i) [format %.8e [lindex $line $i]]
		}

#Schott lists coefficients in the order B1 B2 B3 C1 C2 C3.
#ZEMAX uses B1 C1 B2 C2 B3 C3.  I will use the Schott method, because
#Corning also uses that order.
#Whoops, this works for Sellmeier, but not Cauchy, where coefficients are in
#the standard order.
	   if {$type == "s"} {
		lappend coeffs $coeff(1) $coeff(3) $coeff(5) \
		   $coeff(2) $coeff(4) $coeff(6)
	   } else {
		lappend coeffs $coeff(1) $coeff(2) $coeff(3) \
		   $coeff(4) $coeff(5) $coeff(6)
		}
	   }
	}
   close $fid


#Clean up last filter.
   if {[llength $coeffs] > 0 && [lsearch $names $glass] < 0} {
	puts $outid [format $format $glass $maker $type $coeffs]
	lappend names $glass
	}
   return
   }

#######################################################################
#Read ohara.trans.csv file.
proc oharaDispRead {file outid _innames} {
   upvar 1 $_innames names
   global glassName
   global here
   set format "%-8s %-8s %s %s"
   set maker ohara
   set fid [open $file]
   set nm ""
   set type ""
   set coeffs ""
   gets $fid
   while {1} {
	set line [gets $fid]
	if {[eof $fid]} break
	set line [split $line ,]
	set newname [lindex [lindex $line 0] 0]
	set oldname [lindex [lindex $line 1] 0]

#Sellmeier coeffs
	set a1 [lindex $line 2]
	if {$a1 != ""} {
	   set type s
	   set a2 [lindex $line 3]
	   set a3 [lindex $line 4]
	   set b1 [lindex $line 5]
	   set b2 [lindex $line 6]
	   set b3 [lindex $line 7]
	} else {

#Cauchy coeffs
	   set type c
	   set a0 [lindex $line 8]
	   set a1 [lindex $line 9]
	   set a2 [lindex $line 10]
	   set a3 [lindex $line 11]
	   set a4 [lindex $line 12]
	   set a5 [lindex $line 13]
	   }

#Convert to case sensitive name
	foreach glass "$newname $oldname" {
	   regsub -all " " $glass "" glass
	   if {[info exists glassName($glass)]} {
		set glass $glassName($glass)
	   } else {
		echo No name translation for $glass

#One fixup that I know about.
		regsub LASF $glass LaSF glass
		}
	   if {[lsearch $names $glass] >= 0} continue

	   set coeffs ""

#Schott lists coefficients in the order B1 B2 B3 C1 C2 C3.
#ZEMAX uses B1 C1 B2 C2 B3 C3.  I will use the Schott method, because
#Corning also uses that order.
#Whoops, this works for Sellmeier, but not Cauchy, where coefficients are in
#the standard order.
	   if {$type == "s"} {
		lappend coeffs $a1 $a2 $a3 $b1 $b2 $b3
	   } else {
		lappend coeffs $a0 $a1 $a2 $a3 $a4 $a5
		}
	   puts $outid [format $format $glass $maker $type $coeffs]
	   lappend names $glass
	   }
	}
   close $fid

   return
   }

#######################################################################
#######################################################################
#Read dispersion data
proc glassDispPrep {} {
   global here glassName

#Get translation list
   echo fetching translation list
   opusName
   set outfile $here/DISPERSION.TXT
   set outid [open $outfile w]
   set names ""

#Original Schott - read from ZMAX file.
   set file $here/ZEMAX/SCHOTT.AGF
   zemaxDispRead $file $outid names

#Fake SIO2 - read from ZMAX file.
#Other oddball glasses can be in here as well.

   set file $here/ZEMAX/SIO2.AGF
   zemaxDispRead $file $outid names

#New Schott.  I created a stripped-down .csv from the excel spreadsheet.
   set file $here/ZEMAX/SCHOTT_2000.AGF
   zemaxDispRead $file $outid names

#Old Schott - in case I have old designs
#   set file $here/ZEMAX/OLD_SCHO.AGF
#   zemaxDispRead $file $outid names

#Ohara
   set file $here/ZEMAX/ohara.agf
   zemaxDispRead $file $outid names

#New Ohara
   set file $here/ZEMAX/OHARA_2002.agf
   zemaxDispRead $file $outid names

#Old Ohara - may need for older designs
   set file $here/ZEMAX/OLD_OHAR.AGF
   zemaxDispRead $file $outid names

#Ohara spreadsheet - has still more glasses missed above
   set file $here/../doc/ohara.disp.csv
   oharaDispRead $file $outid names

#Misc - includes F_SILICA
   set file $here/ZEMAX/MISC.AGF
   zemaxDispRead $file $outid names

#Hoya
   set file $here/ZEMAX/HOYA.AGF
   zemaxDispRead $file $outid names

#Sumita
   set file $here/ZEMAX/Sumita.agf
   zemaxDispRead $file $outid names

   close $outid
   return
   }

###########################################################################
###########################################################################

#Density - I have much less information.

########################################################################
#Should have had this before - translate names

proc glassTrans {glass} {
   global glassName
   if {![info exists glassName]} {
	opusName
	}

#Convert to case sensitive name
   if {[info exists glassName($glass)]} {
	set glass $glassName($glass)
   } else {
#	echo No name translation for $glass

#One fixup that I know about.
	regsub LASF $glass LaSF glass
	}
   return $glass
   }

########################################################################
#Read sio2 file.  This is a compressed .csv with just glass name and
#density.
#outid is the output file handle, _innames is a temporary global variable
#to store the ensemble of all glass names.

proc sio2DenRead {file outid _innames} {
   upvar 1 $_innames names
   global here
   set format "%-10s %-10.2f %-10s"
   set maker nature
   set fid [open $file]
   while {1} {
	set line [split [gets $fid] ,]
	if {[eof $fid]} break
	if {[string length $line] == 0} continue
	if {[string index $line 0] == "#"} continue
	set glass [string trim [lindex $line 0] \"]
	set den [lindex $line 1]

#Double-duty names
	set glasses $glass
	if {[regexp {N-(.+)} $glass all oldglass]} {
	   lappend glasses $oldglass
	   }
	foreach glass $glasses {
	   set glass [glassTrans $glass]
	   if {[lsearch $names $glass] >= 0} continue
	   lappend names $glass
	   puts $outid [format $format $glass $den $maker]
	   }
	}
   close $fid
   return
   }

########################################################################
#Read schott file.  This is a compressed .csv with just glass name and
#density.
#outid is the output file handle, _innames is a temporary global variable
#to store the ensemble of all glass names.

proc schottDenRead {file outid _innames} {
   upvar 1 $_innames names
   global here
   set format "%-10s %-10.2f %-10s"
   set maker schott
   set fid [open $file]
   while {1} {
	set line [split [gets $fid] ,]
	if {[eof $fid]} break
	if {[string length $line] == 0} continue
	if {[string index $line 0] == "#"} continue
	set glass [string trim [lindex $line 0] \"]
	set den [lindex $line 1]

#Double-duty names
	set glasses $glass
	if {[regexp {N-(.+)} $glass all oldglass]} {
	   lappend glasses $oldglass
	   }
	foreach glass $glasses {
	   set glass [glassTrans $glass]
	   if {[lsearch $names $glass] >= 0} continue
	   lappend names $glass
	   puts $outid [format $format $glass $den $maker]
	   }
	}
   close $fid
   return
   }

########################################################################
#Read hikari file.  This is a compressed .csv with just glass name and
#density.
#outid is the output file handle, _innames is a temporary global variable
#to store the ensemble of all glass names.

proc hikariDenRead {file outid _innames} {
   upvar 1 $_innames names
   global here
   set format "%-10s %-10.2f %-10s"
   set maker hikari
   set fid [open $file]
   while {1} {
	set line [split [gets $fid] ,]
	if {[eof $fid]} break
	if {[string length $line] == 0} continue
	if {[string index $line 0] == "#"} continue
	set glass [string trim [lindex $line 0] \"]
	set den [lindex $line 1]

#Double-duty names
	set glasses $glass
	if {[regexp {E-(.+)} $glass all oldglass]} {
	   lappend glasses $oldglass
	   }
	foreach glass $glasses {
	   set glass [glassTrans $glass]
	   if {[lsearch $names $glass] >= 0} continue
	   lappend names $glass
	   puts $outid [format $format $glass $den $maker]
	   }
	}
   close $fid
   return
   }

########################################################################
#Read hoya file.  This is a compressed .csv with just glass name and
#density.
#outid is the output file handle, _innames is a temporary global variable
#to store the ensemble of all glass names.

proc hoyaDenRead {file outid _innames} {
   upvar 1 $_innames names
   global here
   set format "%-10s %-10.2f %-10s"
   set maker hikari
   set fid [open $file]
   while {1} {
	set line [split [gets $fid] ,]
	if {[eof $fid]} break
	if {[string length $line] == 0} continue
	if {[string index $line 0] == "#"} continue
	set glass [string trim [lindex $line 0] "\" "]
	set den [lindex $line 1]

#Double-duty names
	set glasses $glass
	if {[regexp {E-(.+)} $glass all oldglass]} {
	   lappend glasses $oldglass
	   }
	foreach glass $glasses {
	   set glass [glassTrans $glass]
	   if {[lsearch $names $glass] >= 0} continue
	   lappend names $glass
	   puts $outid [format $format $glass $den $maker]
	   }
	}
   close $fid
   return
   }

########################################################################
#Read Ohara file.  This is a compressed .csv with just glass name and
#density.
#outid is the output file handle, _innames is a temporary global variable
#to store the ensemble of all glass names.

proc oharaDenRead {file outid _innames} {
   upvar 1 $_innames names
   global here
   set format "%-10s %-10.2f %-10s"
   set maker ohara
   set fid [open $file]
   while {1} {
	set line [split [gets $fid] ,]
	if {[eof $fid]} break
	if {[string length $line] == 0} continue
	if {[string index $line 0] == "#"} continue
	set glass1 [lindex $line 0]
	set glass2 [lindex $line 1]
	set den [lindex $line 2]

#Double-duty names
	set glasses [list $glass1 $glass2]
	foreach glass $glasses {
	   set glass [join [join $glass ""] ""]
	   if {[lsearch $names $glass] >= 0} continue
	   if {[string length $glass] == 0} continue
	   lappend names $glass
	   puts $outid [format $format $glass $den $maker]
	   }
	}
   close $fid
   return
   }

#######################################################################
#Read density data
proc glassDenPrep {} {
   global here glassName


   set outfile $here/DENSITY.TXT
   set outid [open $outfile w]
   set names ""

#Fake SIO2
   set file $here/../doc/sio2-density.csv
   sio2DenRead $file $outid names

#Schott
   set file $here/../doc/schott-density.csv
   schottDenRead $file $outid names

#Ohara
   set file $here/../doc/ohara-density.csv
   oharaDenRead $file $outid names

#Hikari
   set file $here/../doc/hikari-density.csv
   hikariDenRead $file $outid names

#Hoya
   set file $here/../doc/hoya-density.csv
   hoyaDenRead $file $outid names

   close $outid
   return
   }
################################################################
#Fit a Cauchy function to tabulated data.
#This is relatively easy because it is just a polynomial (in the right coords)
#Cauchy coeffs are a0, a1, ... a5
#Let m = 1/wave^2, n2 = n^2
#We have m*n2 = a1 + a0*m + a2*m^2 + ... + a5*m^5

#list is a list of wavelength-index pairs.  Wavelength in microns
proc cauchyFit {list} {
   set xlist ""
   set ylist ""
   foreach pair $list {
	set w [lindex $pair 0]
	set n [lindex $pair 1]
	set m [expr 1./pow($w,2)]
	set n2 [expr pow($n,2)]
	set y [expr $m*$n2]
	lappend xlist $m
	lappend ylist $y
	}
   set coefflist [polyfit $xlist $ylist 6]

#Coefflist from this fit is not in the right order.
   set a1 [lindex $coefflist 0]
   set a0 [lindex $coefflist 1]
   set a2 [lindex $coefflist 2]
   set a3 [lindex $coefflist 3]
   set a4 [lindex $coefflist 4]
   set a5 [lindex $coefflist 5]
   return [list $a0 $a1 $a2 $a3 $a4 $a5]
   }
