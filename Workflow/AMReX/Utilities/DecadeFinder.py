#!/usr/bin/env python3

import numpy as np
import os.path as os
import GlobalVariables.Units    as gvU
import GlobalVariables.Settings as gvS

from Utilities.Files        import GetFileNumberArray, CreatePlotPath
from Utilities.GetFrameData import GetFrameData, GetCentralDensity
from Utilities.Mesh         import CreateLocations

import Utilities.BounceFinder   as BF


global DecadeList
global DecadeTimeList
global DecadeFrameList

DecadeList = 0.0
DecadeTimeList = 0.0
DecadeFrameList = 0


 #=============================================#
#                                               #
#   CreateDecadeData                            #
#                                               #
 #=============================================#
def CreateDecadeData(   PlotDirectories,    \
                        PlotBaseName,       \
                        DataDirectories,    \
                        DataType            ):
    
    global DecadeList
    global DecadeFrameList
    global DecadeTimeList
    
                        
    print( '\n  Creating Density Decade Data' )
    print( '--------------------------------' )
    DecadeList      = [0.0]*gvS.nDirs
    DecadeTimeList  = [0.0]*gvS.nDirs
    DecadeFrameList = [0]*gvS.nDirs
    for i in range(gvS.nDirs):
        Decade, DecadeFrame, DecadeTime = FindDecades(  PlotDirectories[i],    \
                                                        PlotBaseName,          \
                                                        DataDirectories[i],    \
                                                        DataType[i]            )

        DecadeList[i]       = Decade
        DecadeFrameList[i]  = np.array(DecadeFrame)
        DecadeTimeList[i]   = DecadeTime


    return DecadeList, DecadeFrameList, DecadeTimeList
                        
                        
 #=============================================#
#                                               #
#   FindDecades                                 #
#                                               #
 #=============================================#
def FindDecades( PlotDirectory,      \
                 PlotBaseName,       \
                 PathDataDirectory,  \
                 DataType        ):
 
    
    DecadeFile = PathDataDirectory + '/' + '{:}.dat'.format( 'Decades' )
    Check = os.isfile(DecadeFile)

    
    if Check:
    
        print("Reading Decade Information")
        DecadeData = np.loadtxt(DecadeFile)

        Decades     = []
        DecadeFrames= []
        DecadeTimes = []
        
        nDecades = len(DecadeData)
        for i in range(nDecades):
            Decades.append(int(DecadeData[i][0]))
            DecadeFrames.append(int(DecadeData[i][1]))
            DecadeTimes.append(DecadeData[i][2])
        
    
    else:
        print("Calculating Decade Information")
        Decades, DecadeFrames, DecadeTimes = FindDecadeMethod( PlotDirectory,      \
                                                               PlotBaseName,       \
                                                               PathDataDirectory,  \
                                                               DataType,           \
                                                               DecadeFile          )

#    print(Decades)
#    print(DecadeFrames)
#    print(DecadeTimes)
#    exit()
    return Decades, DecadeFrames, DecadeTimes








 #=============================================#
#                                               #
#   FindDecadeMethod                            #
#                                               #
 #=============================================#
def FindDecadeMethod( PlotDirectory,      \
                      PlotBaseName,       \
                      PathDataDirectory,  \
                      DataType,           \
                      DecadeFile          ):

    FileNumberArray = GetFileNumberArray( PlotDirectory,    \
                                          PlotBaseName,     \
                                          -1, -1, 1         )

    NumPltFiles = len(FileNumberArray)

    
    DecadeLimit = 14
    MaxDecade = -1
    Decades = []
    DecadeFrames = []
    DecadeTimes = []
    for i in range(NumPltFiles):

        print( '\r {:}/{:}'.format( i+1, NumPltFiles ), end = '\r' )

        
        PathPlotDirectory = CreatePlotPath( PlotDirectory,      \
                                            PlotBaseName,       \
                                            FileNumberArray[i], \
                                            DataType            )


        Density, DensityUnits, Time = GetCentralDensity( PathPlotDirectory,  \
                                                         DataType            )
            
        CurDecade = int(np.floor(np.log10(Density)))
          

        if CurDecade > MaxDecade:
            MaxDecade = CurDecade
            Decades.append(MaxDecade)
            DecadeFrames.append(FileNumberArray[i])
            DecadeTimes.append(Time.tolist())

        if CurDecade >= DecadeLimit: break
        



    
    
    
    with open( DecadeFile, 'w' ) as FileOut:

        FileOut.write( '# {:}\n'.format( DecadeFile  ) )
        FileOut.write( '# Decade frame_Decade t_Decade[{:}] \n'.format( gvU.TimeUnits),  )

        DecadeLen = len(Decades)
        for i in range(DecadeLen):
            FileOut.write(str(Decades[i])+' '+str(DecadeFrames[i]) + ' ' + str(DecadeTimes[i]) + '\n')


    return Decades, DecadeFrames, DecadeTimes




 #=============================================#
#                                               #
#   FindDecadeMethod                            #
#                                               #
 #=============================================#
def FindDecadeMethod_Bisect( PlotDirectory,      \
                             PlotBaseName,       \
                             PathDataDirectory,  \
                             DataType,           \
                             DecadeFile          ):

    FileNumberArray = GetFileNumberArray( PlotDirectory,    \
                                          PlotBaseName,     \
                                          -1, -1, 1         )

    NumPltFiles = len(FileNumberArray)

    # GetFrameData needs inputs for these values.
    # But, they aren't used when fetching density.
    # So we create theses phonies.
    X1 = np.array([0.0])
    X2 = np.array([0.0])
    X3 = np.array([0.0])
    dX1 = np.array([0.0])
    dX2 = np.array([0.0])
    dX3 = np.array([0.0])

    Decade = 0.0
    for i in range(NumPltFiles):

        print( '\r {:}/{:}'.format( i, NumPltFiles ), end = '\r' )

        PathPlotDirectory = CreatePlotPath( PlotDirectory,      \
                                            PlotBaseName,       \
                                            FileNumberArray[i], \
                                            DataType            )

#        X1, X2, X3, dX1, dX2, dX3, xL, xH               \
#            = CreateLocations(  PathPlotDirectory,      \
#                                TypeIn = "Leaf"      )


        Density, DensityUnits, \
        X1, X2, X3,         \
        dX1, dX2, dX3,      \
        Time                \
            = GetFrameData( PathPlotDirectory,  \
                            DataType,           \
                            'PF_D',             \
                            'True'              )


        CurMaxDensity = max(Density)

        if CurMaxDensity > Decade:
            Decade = CurMaxDensity
            DecadeFrame   = i
            DecadeTime    = Time


    
    with open( DecadeFile, 'w' ) as FileOut:

        FileOut.write( '# {:}\n'.format( DecadeFile  ) )
        FileOut.write( '# frame_Decade t_Decade [{:}] Central Density at Decade [{:}]\n'.format( gvU.TimeUnits, DensityUnits ),  )

        
        FileOut.write( str(DecadeFrame ) + ' ')
        FileOut.write( str(DecadeTime  ) + ' ')
        FileOut.write( str(Decade) + ' ')
        FileOut.write( '\n' )

    return DecadeFrame, DecadeTime, Decade

