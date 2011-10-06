import sys
sys.path.append('/home/jghao/fmccd')
if len(sys.argv) == 1:
    print 'syntax: sub_img imagefile biasfile targetfile'
else:
    from fermiMCCD import *
    
