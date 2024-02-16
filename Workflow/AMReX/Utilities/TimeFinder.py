#!/usr/bin/env python3

import numpy as np
import os.path as os
import GlobalVariables.Units    as gvU
import GlobalVariables.Settings as gvS

from Utilities.Files        import GetFileNumberArray
from Utilities.GetFrameData import GetTimeData



global TimeFrameList

TimeFrameList = 0


 #=============================================#
#                                               #
#   FindTimeFrames                              #
#                                               #
 #=============================================#
def FindTimeFrames( TimeList,           \
                    PlotDirectories,    \
                    PlotBaseName,       \
                    DataType            ):
    


    print( '\n  Finding Time Frames' )
    print( '-----------------------' )
    
    TimeFrameList = [0.0]*gvS.nDirs

    for i in range(gvS.nDirs):
        TimeFrames, TimesFound = FindTimeFrames_1Dir(  TimeList,              \
                                           PlotDirectories[i],    \
                                           PlotBaseName,          \
                                           DataType[i]            )
                                           
        TimeFrameList[i] = TimeFrames
        

    return TimeFrameList


 #=============================================#
#                                               #
#   FindTimeFrames_1Dir                         #
#                                               #
 #=============================================#
def FindTimeFrames_1Dir( TimeList,          \
                         PlotDirectory,     \
                         PlotBaseName,      \
                         DataType           ):


    FileNumberArray = GetFileNumberArray( PlotDirectory,    \
                                          PlotBaseName,     \
                                          -1, -1, 1         )

    NumFiles = len(FileNumberArray)


    if DataType.lower() == 'amrex':
        FirstPlotDirectory = PlotDirectory      \
                          + PlotBaseName        \
                          + '{:}'.format( str(FileNumberArray[0]).zfill(8) )
                          
        LastPlotDirectory = PlotDirectory       \
                          + PlotBaseName        \
                          + '{:}'.format( str(FileNumberArray[-1]).zfill(8) )
    else:
        FirstPlotDirectory = PlotDirectory       \
                          + PlotBaseName[:-4]   \
                          + '_{:}'.format( str(FileNumberArray[0]).zfill(6) )
                          
        LastPlotDirectory = PlotDirectory       \
                          + PlotBaseName[:-4]   \
                          + '_{:}'.format( str(FileNumberArray[-1]).zfill(6) )




    MinTime = GetTimeData( FirstPlotDirectory, DataType )
    MaxTime = GetTimeData( LastPlotDirectory, DataType )
    MinFrameIndex = 0
    MaxFrameIndex = len(FileNumberArray)-1

    NumTimes = len(TimeList)

    TimesFound =[]
    TimeFrame = []
    for t in range(NumTimes):
    
        print( '\r {:}/{:}'.format( t+1, NumTimes ), end = '\r' )
        
        TargetTime = TimeList[t]

        if t > MaxTime:
            print("Time ",t," is outside of the temporal domain.")
            exit()
        elif t < MinTime:
            print("Time ",t," is outside of the temporal domain.")
            exit()

        PrevMinTime = MinTime
        PrevMinIndex = MinFrameIndex
        
        PrevMaxTime = MaxTime
        PrevMaxIndex = MaxFrameIndex


        iter = 0
        Found = False
        while Found == False:
        
            CurIndex = (PrevMinIndex + PrevMaxIndex)//2
            
            CurFrame = FileNumberArray[CurIndex]
            if DataType.lower() == 'amrex':
                CurPlotDirectory = PlotDirectory        \
                                  + PlotBaseName        \
                                  + '{:}'.format( str(CurFrame).zfill(8) )

            else:
                CurPlotDirectory = PlotDirectory        \
                                  + PlotBaseName[:-4]   \
                                  + '_{:}'.format( str(CurFrame).zfill(6) )
                                  
            CurTime = GetTimeData( CurPlotDirectory, DataType )
            
#            print("Iter: ",iter)
#            print(PrevMinIndex,PrevMaxIndex,CurIndex)
#            print(PrevMinTime,PrevMaxTime,CurTime, TargetTime)
#            print("")
            
            
            if TargetTime == CurTime:
                TimeFrame.append(CurFrame)
                TimesFound.append(CurTime)
                Found = True
            elif CurTime in (PrevMinTime,PrevMaxTime):
                TimeFrame.append(CurFrame)
                TimesFound.append(CurTime)
                Found = True
            elif TargetTime < CurTime:
                PrevMaxIndex = CurIndex
                PrevMaxTime  = CurTime
            elif TargetTime > CurTime:
                PrevMinIndex = CurIndex
                PrevMinTime  = CurTime

            iter = iter + 1

    return TimeFrame, TimesFound




