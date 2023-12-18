#!/usr/bin/env python3

import GlobalVariables.Settings as gvS


 #=============================================#
#                                               #
#   CoordSystemCheck                            #
#                                               #
 #=============================================#
def CoordSystemCheck():

    msg = 'Invalid choice of CoordinateSystem: {:s}'.format( gvS.CoordinateSystem )
    msg += '\n\nValid Choices:\n'
    msg +=   '--------------\n'
    msg +=   'cartesian\n'
    msg +=   'cylindrical\n'
    msg +=   'spherical'

    assert (    gvS.CoordinateSystem == 'cartesian' \
             or gvS.CoordinateSystem == 'cylindrical' \
             or gvS.CoordinateSystem == 'spherical' ), msg

    return
