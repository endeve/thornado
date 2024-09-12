#!/usr/bin/env python3

import numpy as np
import os.path as os
import GlobalVariables.Units    as gvU
import GlobalVariables.Settings as gvS

from Utilities.Files        import GetFileNumberArray
from Utilities.GetFrameData import GetFrameData, GetCentralDensity
from Utilities.Mesh         import CreateLocations



global BounceDensityList
global BounceTimeList
global BounceFrameList

BounceDensityList = 0.0
BounceTimeList = 0.0
BounceFrameList = 0


 #=============================================#
#                                               #
#   CreateBounceData                            #
#                                               #
 #=============================================#
def CreateBounceData(   PlotDirectories,    \
                        PlotBaseName,       \
                        DataDirectories,    \
                        DataType            ):
    
    global BounceDensityList
    global BounceTimeList
    global BounceFrameList
    
                        
    print( '\n  Creating Bounce Data' )
    print( '------------------------' )
    BounceDensityList = [0.0]*gvS.nDirs
    BounceTimeList    = [0.0]*gvS.nDirs
    BounceFrameList   = [0]*gvS.nDirs
    for i in range(gvS.nDirs):
        BFrame, BTime, BDensity = FindBounce(   PlotDirectories[i],    \
                                                PlotBaseName,          \
                                                DataDirectories[i],    \
                                                DataType[i]            )

        BounceDensityList[i] = BDensity
        BounceFrameList[i]   = BFrame
        BounceTimeList[i]    = BTime


                        
                        
                        
 #=============================================#
#                                               #
#   FindBounce                                  #
#                                               #
 #=============================================#
def FindBounce( PlotDirectory,      \
                PlotBaseName,       \
                PathDataDirectory,  \
                DataType        ):


    BounceFile = PathDataDirectory + '/' + '{:}.dat'.format( 'Bounce' )
    Check = os.isfile(BounceFile)
    
    if Check:
    
        print("Reading Bounce Information")
        Bounce = np.loadtxt(BounceFile)
        
        
        BounceFrame = int(Bounce[0])
        BounceTime = Bounce[1]
        BounceDensity = Bounce[2]
    
    else:
        print("Calculating Bounce Information")
        BounceFrame, BounceTime, BounceDensity = FindBounceMethod( PlotDirectory,      \
                                                                   PlotBaseName,       \
                                                                   PathDataDirectory,  \
                                                                   DataType,           \
                                                                   BounceFile          )


    return BounceFrame, BounceTime, BounceDensity








 #=============================================#
#                                               #
#   FindBounceMethod                            #
#                                               #
 #=============================================#
def FindBounceMethod( PlotDirectory,      \
                      PlotBaseName,       \
                      PathDataDirectory,  \
                      DataType,           \
                      BounceFile          ):

    FileNumberArray = GetFileNumberArray( PlotDirectory,    \
                                          PlotBaseName,     \
                                          -1, -1, 1         )

    NumPltFiles = len(FileNumberArray)
    


    BounceFlagTol = 100
    BounceFlag = 0
    BounceDensity = 0
    BounceDensityThresholdFlag = 0
    BounceDensityThreshold = 1.0E12
    BounceFrame = -1
    BounceTime = -1.0
    for i in range(NumPltFiles):

        print( '\r {:}/{:}'.format( i, NumPltFiles ), end = '\r' )

        PlotFileNumber    = FileNumberArray[i]
        if DataType.lower() == 'amrex':
            PathPlotDirectory = PlotDirectory       \
                              + PlotBaseName    \
                              + '{:}'.format( str(PlotFileNumber).zfill(8) )
        else:
            PathPlotDirectory = PlotDirectory       \
                              + PlotBaseName[:-4]   \
                              + '_{:}'.format( str(PlotFileNumber).zfill(6) )

                            
        Density, DensityUnits, Time = GetCentralDensity( PathPlotDirectory,  \
                                                         DataType            )
                                                         
          
#        CurMaxDensity = max(Density)
        CurMaxDensity = Density[0]
        
#        print(CurMaxDensity,BounceDensity,Time,i,BounceFlag,BounceDensityThresholdFlag)

        if not BounceDensityThresholdFlag:
            if CurMaxDensity > BounceDensityThreshold:
               BounceDensityThresholdFlag = True
               
        if BounceDensityThresholdFlag:
            if CurMaxDensity > BounceDensity:
                BounceDensity = CurMaxDensity
                BounceFrame   = i
                BounceTime    = Time
                BounceFlag    = 0
            else:
                BounceFlag = BounceFlag + 1
                if BounceFlag > BounceFlagTol: break
      
  
    if BounceFrame == -1:
        print('Bounce not found.  Exiting')
        exit()
    
    
    else:
        with open( BounceFile, 'w' ) as FileOut:

            FileOut.write( '# {:}\n'.format( BounceFile  ) )
            FileOut.write( '# frame_bounce t_Bounce [{:}] Central Density at bounce [{:}]\n'.format( gvU.TimeUnits, DensityUnits ),  )

            
            FileOut.write( str(BounceFrame ) + ' ')
            FileOut.write( str(BounceTime  ) + ' ')
            FileOut.write( str(BounceDensity) + ' ')
            FileOut.write( '\n' )

    return BounceFrame, BounceTime, BounceDensity




 #=============================================#
#                                               #
#   FindBounceMethod                            #
#                                               #
 #=============================================#
def FindBounceMethod_Bisect( PlotDirectory,      \
                             PlotBaseName,       \
                             PathDataDirectory,  \
                             DataType,           \
                             BounceFile          ):

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

    BounceDensity = 0.0
    for i in range(NumPltFiles):

        print( '\r {:}/{:}'.format( i, NumPltFiles ), end = '\r' )

        PlotFileNumber    = FileNumberArray[i]
        if DataType.lower() == 'amrex':
            PathPlotDirectory = PlotDirectory       \
                              + PlotBaseName    \
                              + '{:}'.format( str(PlotFileNumber).zfill(8) )
        else:
            PathPlotDirectory = PlotDirectory       \
                              + PlotBaseName[:-4]   \
                              + '_{:}'.format( str(PlotFileNumber).zfill(6) )

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

        if CurMaxDensity > BounceDensity:
            BounceDensity = CurMaxDensity
            BounceFrame   = i
            BounceTime    = Time


    
    with open( BounceFile, 'w' ) as FileOut:

        FileOut.write( '# {:}\n'.format( BounceFile  ) )
        FileOut.write( '# frame_bounce t_Bounce [{:}] Central Density at bounce [{:}]\n'.format( gvU.TimeUnits, DensityUnits ),  )

        
        FileOut.write( str(BounceFrame ) + ' ')
        FileOut.write( str(BounceTime  ) + ' ')
        FileOut.write( str(BounceDensity) + ' ')
        FileOut.write( '\n' )

    return BounceFrame, BounceTime, BounceDensity
