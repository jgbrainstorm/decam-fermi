#Put a string at x,y then return cursor to lower left corner
proc vtPutxy {x y char} {
   puts stdout \033\[${x}\;${y}H${char} nonewline
#   puts stdout \033\[24\;0H nonewline
   }

proc vtMovexy {x y} {
   puts stdout \033\[${x}\;${y}H nonewline
   }

#Set scroll region from n1 to n2
proc vtScroll {n1 n2} {
   puts stdout \033\[${n1}\;${n2}r nonewline
   puts stdout \033\[${n2}\;0H nonewline
   return
   }

proc vtBold {} {
   puts stdout \033\[1m nonewline
   return
   }

proc vtReverse {} {
   puts stdout \033\[7m nonewline
   return
   }

proc vtNormal {} {
   puts stdout \033\[m nonewline
   return
   }

proc vtSaveCursor {} {
   puts stdout \0337 nonewline
   return
   }

proc vtRestoreCursor {} {
   puts stdout \0338 nonewline
   return
   }

proc vtInsertLine {} {
   puts stdout \033\[L nonewline
   return
   }

proc vtClearToDisplayEnd {} {
   puts stdout \033\[J nonewline
   return
   }

proc vtClearToLineEnd {} {
   puts stdout \033\[K nonewline
   return
   }

proc vtClearScreen {} {
   puts stdout \033\[H\033\[2J nonewline
   return
   }

proc screen {} {
   vtClearScreen
   vtMovexy 24 1
   return
   }

proc vtDelChar {{n 1}} {
   puts stdout \033\[${n}P nonewline
   return
   }

proc vtDelLine {} {
   puts stdout \033\[M nonewline
   return
   }

proc vtHome {} {
   puts stdout \033\[H nonewline
   return
   }

proc vtReset {} {
# :is=\E[r\E[m\E[2J\E[H\E[?7h\E[?1;3;4;6l\E[4l:\
   puts stdout \033\[r nonewline
   puts stdout \033\[m nonewline
   puts stdout \033\[2J nonewline
   puts stdout \033\[H nonewline
   puts stdout \033\[7h nonewline
   puts stdout \033\[1 nonewline
   puts stdout \033\[3 nonewline
   puts stdout \033\[4 nonewline
   puts stdout \033\[5 nonewline
   puts stdout \033\[4l
   return
   }

proc vtUnderline {} {
   puts stdout \033\[4m nonewline
   return
   }

proc vtBell {} {
   puts stdout [format %c 7] nonewline 
   return
   }

#The following is a template for printing lots of vt100 special symbols
#char can be a-z and some punctuation marks.  Unused chars are echoed normally.

proc vtBox {} {
   set char a
   puts stdout [format \033\[0m%c${char}\033\[0m%c 14 15] nonewline
   }

proc vtAlt {line {nonewline ""}} {
   if {$nonewline == ""} then {
	puts stdout [format \033\[0m%c${line}\033\[0m%c 14 15]
   } else {
	puts stdout [format \033\[0m%c${line}\033\[0m%c 14 15] nonewline
	}
   }

#Use normal char set
proc vtSet0 {} {
   puts stdout [format \033\[0m%c 15] nonewline
   }

#Use alternate char set
proc vtSet1 {} {
   puts stdout [format \033\[0m%c 14] nonewline
   }

#Move cursor around on screen by n chars
proc vtLeft {{n 1}} {
   if {$n == 0} return
   puts stdout \033\[${n}D nonewline
   }

proc vtRight {{n 1}} {
   if {$n == 0} return
   puts stdout \033\[${n}C nonewline
   }

proc vtUp {{n 1}} {
   if {$n == 0} return
   puts stdout \033\[${n}A nonewline
   }

proc vtDown {{n 1}} {
   if {$n == 0} return
   puts stdout \033\[${n}B nonewline
   }

proc vtWrap {n} {
   if {$n == 0} {
	puts stdout \033\[?7l nonewline
	puts stdout \033\[?45l nonewline
	}
   if {$n == 1} {
	puts stdout \033\[?7h nonewline
	puts stdout \033\[?45h nonewline
	}
   }

proc vtWrapSave {} {
   puts stdout \033\[?7s nonewline
   puts stdout \033\[?45s nonewline
   }

proc vtWrapRestore {} {
   puts stdout \033\[?7r nonewline
   puts stdout \033\[?45r nonewline
   }

proc vtBlank {{n 1}} {
   puts stdout \033\[${n}@ nonewline
   }

