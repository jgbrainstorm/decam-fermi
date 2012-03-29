#Load auxiliary files from same directory as ile.tcl
set script [info script]
set iledir [file dirname $script]

if {[info command loop] == ""} {source $iledir/tclx.tcl}
source $iledir/vt100.tcl
unset iledir

if {[file exists $tcl_library/parray.tcl]} {source $tcl_library/parray.tcl}

#Move cursor left or right, including bounds checking. Input variable ichar
#is passed by reference and updated
proc curMove {nmove line inichar} {
   upvar 1 $inichar ichar
   set len [string length $line]
   if {$nmove < 0} {
	set nmove [max $nmove [expr -1*$ichar]]
	vtLeft [expr -1*$nmove]
	}
   if {$nmove > 0} {
	set nmove [min $nmove [expr $len-$ichar]]
	set chars [string range $line $ichar [expr $ichar+$nmove-1]]
	puts stdout $chars nonewline
	}
   incr ichar $nmove
   return
   }

#Fetch the next character.  We allow an external line to be input.  This can
#be used for a type-ahead buffer.

proc getNextChar {} {
   global TYPEAHEAD
   if {[info exists TYPEAHEAD]} {
	if {[string length $TYPEAHEAD] > 0} {
	   set char [string index $TYPEAHEAD 0]
	   set TYPEAHEAD [string range $TYPEAHEAD 1 end]
	   return $char
	   }
	}
#Create a little polling loop to read the next character while permitting update
#operations.  There is probably a more elegant way to do this.

   while {1} {

#Use the select system call to test for new input; also wait .1 seconds
#If I use the tcl "after" command, it causes a freeze after a laptop suspend.
	set ready [select stdin "" "" .1]
	if {$ready != ""} {
	   set char [read stdin 1]
	   if {[string length $char] > 0} break
	   }
#Carry out update of screen and other local events (if command exists)
	if {[info command update] != ""} update
	if {[info command localupdate] != ""} localupdate
	}
   return $char
   }

#Fetch next character and store in TYPEAHEAD buffer.  Return 0 if normal
#character, 1 if a control-C

proc getTypeAhead {} {
   global TYPEAHEAD
#Get current blocking mode
   set block [fconfigure stdin -blocking]
   fconfigure stdin -blocking 0
   set return 0
   while {1} {
	set char [read stdin 1]
	if {[fblocked stdin]} break
#Check if character is control-C or -D
	scan $char %c dec
	if {$dec == 3 || $dec == 4} {
	   if {[info exists TYPEAHEAD]} {unset TYPEAHEAD}
	   set return $dec
	   break
	   }
	append TYPEAHEAD $char
	}
   fconfigure stdin -blocking $block
   return $return
   }

#Read a line
proc getLine {{prompt "prompt> "}} {
#   globals
#Set autowrap and reverse autowrap
   vtWrapSave
   vtWrap 1
   set line ""
   exec stty raw -echo
#Which line in history record?
   set ihist [history nextid]
   set nhist $ihist
#Current history line is stored separately
   set curhistline ""
   set ESC [format %c 033]
   set ichar 0
   set ndel 0
   set erase ""
#Print out the prompt
   puts stdout $prompt nonewline
#Prepare for multiline input
   set startline ""
   while {1} {
#ichar is current location of cursor in line. imove is the number of spaces to
#move the cursor from the previously typed character.  ndel is the number of
#characters deleted; we will append spaces to erase any residual.
#New cursor position
#If erase > 0, clear out old line and redisplay new line
	if {$erase != ""} {
	   vtLeft $ichar
	   for {set i 0} {$i < $erase} {incr i} {puts stdout " " nonewline}
	   vtLeft $erase
	   puts stdout $line nonewline
	   set ichar [string length $line]
	   }
	set len [string length $line]
#Cursor should be over 1st char of part2.  Reparse line into 2 parts
	set part1 [string range $line 0 [expr $ichar-1]]
	set part2 [string range $line $ichar end]
#Redisplay part2
	puts stdout $part2 nonewline
	for {set i 0} {$i < $ndel} {incr i} {
	   puts stdout " " nonewline
	   }
	vtLeft [expr [string length $part2]+$ndel]
	flush stdout
#Reset some params
	set ndel 0
	set erase ""
#After processing any character, the following state variables must be correct:
#   ichar is the current cursor position
#   line is the full text of the line
#   ndel is the number of characters at the end of the line to be replaced with
#	blanks
#   erase is the number of characters in a line (from the beginning) to be replaced
#	with blanks.

#Full line is startline plus current line
	set fullline $startline$line

#Now read a new character
	set char [getNextChar]

#Delete (BS or Del)
	if {$char == "\b" || $char == [format %c 127]} {
	   if {$ichar == 0} continue
	   set part1  [string range $part1 0 [expr [string length $part1]-2]]
	   set ndel 1
	   curMove -1 $line ichar
#Return or newline, which seems to be the tcl7.5 version of <cr>
	} elseif {$char == "\r" || $char == "\n"} {
	   if {[info complete $fullline]} {
		break
	   } else {
		set startline $fullline\n
#puts stdout $curhistline
		if {$line != ""} {history add $line}
		set ihist [history nextid]
		set nhist $ihist
		set line ""
		puts stdout "\r\n===>" nonewline
		set part1 ""; set part2 ""
		set ichar 0
		}
	} elseif {$char == $ESC} then {
           set char2 [read stdin 1]
           set char3 [read stdin 1]
           if {$char2 != "\["} then continue
#Up arrow
           if {$char3 == "A"} then {
		if {[catch {history event [expr $ihist-1]}]} continue
		incr ihist -1
		set erase [string length $line]
		set part1 [history event $ihist]
		set part2 ""
                }
#Down arrow
           if {$char3 == "B"} then {
		if {$ihist < $nhist} {incr ihist}
		set erase [string length $line]
		if {$ihist < $nhist} {
		   set part1 [history even $ihist]
		} else {
		   set part1 $curhistline
		   }
		set part2 ""
                }
#Right arrow
           if {$char3 == "C"} then {
		curMove 1 $line ichar
                }
#Left arrow
         if {$char3 == "D"} then {
                curMove -1 $line ichar
                }
#Ctrl-E
	} elseif {$char == "[format %c 5]"} {
	   curMove [expr [string length $line] - $ichar] $line ichar
#Ctrl-K
	} elseif {$char == "[format %c 11]"} {
	   set ndel [string length $part2]
	   set part2 ""
#Ctrl-U
	} elseif {$char == "[format %c 21]"} {
	   curMove [expr -1*[string length $part1]] $line ichar
	   set ndel [string length $part1]
	   set part1 ""
#Printable character
	} elseif {[ctype print $char]} {
	   append part1 $char
	   puts stdout $char nonewline
	   flush stdout
	   incr ichar
#Unknown
	} else {
	   break
	   }
	set line $part1$part2
	if {$ihist < $nhist} {
#	   catch {history change $line $ihist}
	   continue
	   }
	if {$ihist == $nhist} {
	   set curhistline $line
	   }
	}
   set curhistline $line
   if {$curhistline != ""} {history add $curhistline}
   set line $fullline
   exec stty cooked echo
   puts stdout ""
   flush stdout
   vtWrapRestore
   return $line
   }

#Put out a line to the tty when in raw mode

proc putLine {text} {
   puts stdout \r${text}\r
   return
   }

proc ile {} {
   global tcl_prompt1
   while {1} {
	set prog [info nameofexecutable]
	set prompt "[file tail $prog]> "
	set line [getLine $prompt]
	if {[string index $line 0] == "!"} {
	   set pattern [string range $line 1 end]
	   if {[catch {set line [history event $pattern]}]} {
		puts stdout "Error: No match for $pattern"
		continue
	   	}
	   global TYPEAHEAD
	   set TYPEAHEAD $line
	   continue
	   }
#	if {[catch {uplevel #0 $line} msg]} {catch {eval unknown $line} msg}
#Extract the command name and see if it is internal.  Don't assume that "line" is
#a valid TCL list.
	set list [split $line]
	set command [string trimright [lindex $list 0] \;]
	if {$command == ""} continue
	if {[info command $command] != ""} {
	   catch {uplevel #0 $line} msg
	   if {$msg != ""} {putLine $msg}
	} else {
            set check [auto_execok $command]
            if {$check != ""} {
		catch {uplevel #0 exec sh -c \"$line\" >&@ stdout} msg
		if {$msg != ""} {putLine $msg}
	   } else {
		echo Unknown command: $line
		}
	   }
	}
   }

#If we are not in interactive mode, don't run all the stuff below.
#In particular, I will not try to source yet another file based on argv,
#since in non-interactive mode I like to use argv to pass arguments to other
#procs.

if {![fstat stdin tty]} return

#Withdraw that stupid . window
if {[info command .] != ""} {
	wm withdraw .
	}

set dir [file dir [info script]]

if {![info complete $argv]} {
   echo Command line error: $argv
   exit
   }

if {[string compare [lindex $argv 0] "--"] == 0} {
   set argv [lrange $argv 1 end]
} else {
   if {$argv != ""} {
	set cmdfile [lindex $argv 0]
	if {[catch {source $cmdfile} msg]} {puts stdout $msg}
	}
   }

flush stdout

ile
