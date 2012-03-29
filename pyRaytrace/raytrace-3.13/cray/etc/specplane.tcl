proc specflag {hndl} {
   setparam $hndl 10 1 0
   setparam $hndl 10 8 0
   setparam $hndl 10 9 0
   setparam $hndl 10 10 0
   setparam $hndl 10 11 0
   setflag $hndl 10 5 1
   setflag $hndl 10 8 1
   setflag $hndl 10 9 1
   setflag $hndl 10 10 1
   setflag $hndl 0 2 2
   setflag $hndl 0 3 2
   setflag $hndl 0 4 2
   setflag $hndl 0 5 2
   }
