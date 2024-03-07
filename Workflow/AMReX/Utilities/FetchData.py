#!/usr/bin/env python3

import numpy as np
from sys import argv



 #=============================================#
#                                               #
#   fetchData_AMReX                             #
#                                               #
 #=============================================#
def fetchData_AMReX(t, FileNumberArray, DataDirectory, Field ):

    FileDirectory = DataDirectory + str(FileNumberArray[t]) + '/'

    TimeFile = FileDirectory + '{:}.dat'.format( 'Time' )
    X1File   = FileDirectory + '{:}.dat'.format( 'X1' )
    dX1File  = FileDirectory + '{:}.dat'.format( 'dX1' )
    DataFile = FileDirectory + '{:}.dat'.format( Field )

    DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )

    Time = np.loadtxt( TimeFile )
    X1_C = np.loadtxt( X1File   )
    dX1  = np.loadtxt( dX1File  )
    Data = np.loadtxt( DataFile )

    return Data, DataUnits, X1_C, dX1, Time


 #=============================================#
#                                               #
#   fetchData_Native                             #
#                                               #
 #=============================================#
def fetchData_Native(t, FileNumberArray, DataDirectory, Field ):


    FileDirectory = DataDirectory + str(FileNumberArray[t]) + '/'

    TimeFile = FileDirectory + '{:}.dat'.format( 'Time' )
    X1File   = FileDirectory + '{:}.dat'.format( 'X1' )
    dX1File  = FileDirectory + '{:}.dat'.format( 'dX1' )
    DataFile = FileDirectory + '{:}.dat'.format( Field )

    DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )

    Time = np.loadtxt( TimeFile )
    X1_C = np.loadtxt( X1File   )
    dX1  = np.loadtxt( dX1File  )
    Data = np.loadtxt( DataFile )

    return Data, DataUnits, X1_C, dX1, Time




 #=============================================#
#                                               #
#   ReadHeader                                  #
#                                               #
 #=============================================#
def ReadHeader( DataFile ):

    f = open( DataFile )

    dum = f.readline()

    s = f.readline(); ind = s.find( ':' )+1
    DataShape = np.array( list( map( np.int64, s[ind:].split() ) ), np.int64 )

    s = f.readline(); ind = s.find( ':' )+1
    DataUnits = s[ind:]

    s = f.readline(); ind = s.find( ':' )+1
    MinVal = np.float64( s[ind:] )

    s = f.readline(); ind = s.find( ':' )+1
    MaxVal = np.float64( s[ind:] )

    f.close()

    return DataShape, DataUnits, MinVal, MaxVal




 #=============================================#
#                                               #
#   fetchTime_AMReX                             #
#                                               #
 #=============================================#
def fetchTime_AMReX(t, FileNumberArray, DataDirectory ):

    FileDirectory = DataDirectory + str(FileNumberArray[t]) + '/'

    TimeFile = FileDirectory + '{:}.dat'.format( 'Time' )

    Time = np.loadtxt( TimeFile )

    return Time



 #=============================================#
#                                               #
#   fetchData_AMReX                             #
#                                               #
 #=============================================#
def fetchCenDen_AMReX(t, FileNumberArray, DataDirectory ):

    FileDirectory = DataDirectory + str(FileNumberArray[t]) + '/'

    DataFile = FileDirectory + '{:}.dat'.format( 'PF_D' )

    DataShape, DataUnits, MinVal, MaxVal = ReadHeader( DataFile )
    

    return MaxVal
