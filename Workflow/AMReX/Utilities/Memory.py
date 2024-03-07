#!/usr/bin/env python3

#=============================================#
#   Included Routines
#
#   PrintMemory
#
#=============================================#


 #=============================================#
#                                               #
#   PrintMemory                                 #
#                                               #
 #=============================================#
def PrintMemory()
    import os
    import psutil
    process = psutil.Process( os.getpid() )
    print( 'mem: {:.3e} kB'.format \
            ( process.memory_info().rss / 1024.0 ) )
