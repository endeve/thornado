#!/usr/bin/env python3

import numpy as np
from sys import argv

import GlobalVariables.Settings as gvS
import Utilities.BounceFinder   as BF

from Utilities.FetchData import fetchTime_AMReX, fetchCenDen_AMReX

#=============================================#
#   Included Routines
#
#   SychronizeFrameList
#
#=============================================#
def SynchronizeFrameLists(nDirs,nSS, FileNumberArray, DataDirectory, SyncBy):
        
    if SyncBy == 'Time':
        TimeList = [];
        for i in range(nDirs):
            Time_List.append(['None']*nSS[i])
            for j in range(nSS[i]):
                
                Time = fetchTime_AMReX( j,                  \
                                        FileNumberArray[i], \
                                        DataDirectory[i]    )
                TimeList[i][j] = Time
        
        
        SubFrameList = Synchronizer_Time(nDirs,TimeList,nSS)
        
    elif SyncBy == 'Bounce':
    
        # Check for bounce data
        if ( BF.BounceTimeList == 0.0 ):
            exit("Must be run with 'gvS.ReferenceBounce = True'. ")
    


        DensityList = [];
        for i in range(nDirs):
            DensityList.append(['None']*BF.BounceFrameList[i])
            for j in range(BF.BounceFrameList[i]):
                Density = fetchCenDen_AMReX( j,                  \
                                             FileNumberArray[i], \
                                             DataDirectory[i]    )
                DensityList[i][j] = Density

        SubFrameListA = Synchronizer_Density(nDirs,DensityList,BF.BounceFrameList)
        SubFrameLengthA = len(SubFrameListA)



        nSS_Time = [0.0]*nDirs
        TimeList = []
        for i in range(nDirs):
            nSS_Time[i] = nSS[i]-BF.BounceFrameList[i]
            TimeList.append(['None']*nSS_Time[i])
            for j in range(nSS_Time[i]):
                Time = fetchTime_AMReX(BF.BounceFrameList[i]+j, \
                                       FileNumberArray[i], \
                                       DataDirectory[i]    )
                
                TimeList[i][j] = Time
                
        SubFrameListB = Synchronizer_Time(nDirs,TimeList,nSS_Time)
        SubFrameLengthB = len(SubFrameListB)

        for i in range(nDirs):
            for j in range(SubFrameLengthB):
                SubFrameListB[j][i] = SubFrameListB[j][i] + BF.BounceFrameList[i]

        SubFrameList = [];
        for i in range(SubFrameLengthA):
            SubFrameList.append(SubFrameListA[i][:])
        for i in range(SubFrameLengthB):
            SubFrameList.append(SubFrameListB[i][:])
        
        
    SubFrameLength = len(SubFrameList)

    FrameList = [];
    PrevFrame = [0]*nDirs
    FrameList.append(PrevFrame)
    for i in range(SubFrameLength):
        FrameList.append(SubFrameList[i][:])




#    print('FrameList')
#    for i in range(len(FrameList)):
#        print(FrameList[i])
#    exit()
    return FrameList


 #=============================================#
#                                               #
#   Synchronizer_Time                           #
#                                               #
 #=============================================#
def Synchronizer_Time(nDirs,TimeLists,nSS):

    SyncTol = 1.0E-3

    FrameList = [];
    PrevFrame = [0]*nDirs
    
#   Sort the lists until you reach the end of one.
    ListsFinished = 0
    while ListsFinished == 0:
    
#       Find the current minimum time
        CurTimes = []
        for i in range(nDirs):
            CurTimes.append(TimeLists[i][PrevFrame[i]])
        MinTime = min(CurTimes)


#       Identify all times within a tolerance of the min.
        MinTimeListIndex = []
        for i in range(nDirs):
            Time = TimeLists[i][PrevFrame[i]]
            if abs(Time - MinTime) <= SyncTol*MinTime:
                MinTimeListIndex.append(i)
        
        
#       Advance the frames with minimum time
        for i in range(len(MinTimeListIndex)):
            j = MinTimeListIndex[i]
            PrevFrame[j] = min(PrevFrame[j]+1,nSS[j])
            

#       Add frames to list
        FrameList.append(PrevFrame[:])
        
        
#       Detect if a list is at it's end.
        EoL_Flag = [0]*nDirs
        for i in range(nDirs):
            if PrevFrame[i] == nSS[i]-1:
                EoL_Flag[i] = 1
                
                
#       Detect if a list has reached the stop time.
        for i in range(nDirs):
            if EoL_Flag[i] == 0:
                curTime = TimeLists[i][PrevFrame[i]]
                if gvS.ReferenceBounce:
                    curTime = curTime-BF.BounceTimeList[i]
                if curTime > gvS.StopTime:
                    EoL_Flag[i] = 1

#       If ListsFinished > 0, while loop ends.
        ListsFinished = sum(EoL_Flag)
        
 
#   Create new lists with finished list(s) removed and start recurusion.
    LastFrame = PrevFrame[:]
    if ListsFinished != nDirs:
        NewTimeLists = []
        NewnSS = []
        OldToNewMap = [-1]*nDirs
        NewIndex = 0
        for i in range(nDirs):
            if EoL_Flag[i] == 0:
                NewTimeLists.append(TimeLists[i][PrevFrame[i]:nSS[i]])
                NewnSS.append(nSS[i]-(PrevFrame[i]))
                OldToNewMap[i] = NewIndex
                NewIndex += 1
                
        NewnDirs = nDirs-ListsFinished
        
        
#       Recursive call.
        SubFrameList = Synchronizer_Time(NewnDirs,NewTimeLists,NewnSS)

 
#       Add returned list to FrameList and fill in frames for finished list(s).
        SubFrameLength = len(SubFrameList)
        for i in range(SubFrameLength):
            for j in range(nDirs):
                if EoL_Flag[j] == 0:
                    PrevFrame[j] = SubFrameList[i][OldToNewMap[j]]+LastFrame[j]
                elif EoL_Flag[j] == 1:
                    PrevFrame[j] = PrevFrame[j]
                
            FrameList.append(PrevFrame[:])

        
    return FrameList









 #=============================================#
#                                               #
#   Synchronizer_Density                        #
#                                               #
 #=============================================#
def Synchronizer_Density(nDirs,DensityLists,nSS):

    SyncTol = 1.0E-3

    FrameList = [];
    PrevFrame = [0]*nDirs

    
#   Sort the lists until you reach the end of one.
    ListsFinished = 0
    while ListsFinished == 0:

#       Find the current minimum density
        CurDensity = []
        for i in range(nDirs):
            CurDensity.append(DensityLists[i][PrevFrame[i]])
        MinDen = min(CurDensity)
        
        
        
#       Identify all densities within a tolerance of the min.
        MinDenListIndex = []
        for i in range(nDirs):
            
            Density = DensityLists[i][PrevFrame[i]]
            if abs(Density - MinDen) <= SyncTol*MinDen:
                MinDenListIndex.append(i)
        
        
        
        
#       Advance the frames with minimum density
        for i in range(len(MinDenListIndex)):
            j = MinDenListIndex[i]
            PrevFrame[j] = min(PrevFrame[j]+1,nSS[j])
            

#       Add frames to list
        FrameList.append(PrevFrame[:])
        

#       Detect if a list is at it's end.
        EoL_Flag = [0]*nDirs
        for i in range(nDirs):
            if PrevFrame[i] == nSS[i]-1:
                EoL_Flag[i] = 1
                

#       If ListsFinished > 0, while loop ends.
        ListsFinished = sum(EoL_Flag)
        
 
#   Create new lists with finished list(s) removed and start recurusion.
    LastFrame = PrevFrame[:]
    if ListsFinished != nDirs:
        NewDenLists = []
        NewnSS = []
        OldToNewMap = [-1]*nDirs
        NewIndex = 0
        for i in range(nDirs):
            if EoL_Flag[i] == 0:
                NewDenLists.append(DensityLists[i][PrevFrame[i]:nSS[i]])
                NewnSS.append(nSS[i]-(PrevFrame[i]))
                OldToNewMap[i] = NewIndex
                NewIndex += 1
                
        NewnDirs = nDirs-ListsFinished
        
        
#       Recursive call.
        SubFrameList = Synchronizer_Density(NewnDirs,NewDenLists,NewnSS)

 
#       Add returned list to FrameList and fill in frames for finished list(s).
        SubFrameLength = len(SubFrameList)
        for i in range(SubFrameLength):
            for j in range(nDirs):
                if EoL_Flag[j] == 0:
                    PrevFrame[j] = SubFrameList[i][OldToNewMap[j]]+LastFrame[j]
                elif EoL_Flag[j] == 1:
                    PrevFrame[j] = PrevFrame[j]
                
            FrameList.append(PrevFrame[:])

        
    return FrameList
