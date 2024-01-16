#!/usr/bin/env python3


global TimeUnits
global X1Units
global X2Units
global X3Units

# Default Units 
TimeUnits = ''
X1Units   = ''
X2Units   = ''
X3Units   = ''



def SetSpaceTimeUnits(CoordinateSystem, UsePhysicalUnits):

    global TimeUnits
    global X1Units
    global X2Units
    global X3Units
    TimeUnits = ''
    X1Units   = ''
    X2Units   = ''
    X3Units   = ''
    if UsePhysicalUnits:
        TimeUnits = 'ms'
        if   CoordinateSystem == 'cartesian':
            X1Units = 'km'
            X2Units = 'km'
            X3Units = 'km'
        elif CoordinateSystem == 'spherical':
            X1Units = 'km'
            X2Units = 'rad'
            X3Units = 'rad'
            
    return
