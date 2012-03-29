#
#Initialization file
#
global env
set dir [file dirname [info script]]
source $dir/shTools.tcl
source $dir/shGlobal.tcl
source $dir/shParse.tcl
source $dir/shHdr.tcl
source $dir/shHandle.tcl
source $dir/shSchema.tcl
source $dir/shP2C.tcl
source $dir/screenEdit.tcl
source $dir/shCurses.tcl
source $dir/vt100.tcl
source $dir/shChain.tcl
source $dir/pg.tcl
source $dir/shMemory.tcl

#No other place to declare this - increase the size of the history list.
history keep 1000
return
