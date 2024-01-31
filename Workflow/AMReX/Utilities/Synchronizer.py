#!/usr/bin/env python3

import numpy as np
from sys import argv

import GlobalVariables.Settings as gvS
import Utilities.BounceFinder   as BF

#=============================================#
#   Included Routines
#
#   SychronizeFrameList
#
#=============================================#
def SynchronizeFrameLists(nDirs,nSS,SyncBy):

    
    
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
    else if SyncBy == 'Density':
        DensityList = [];
        for i in range(nDirs):
            Density_List.append(['None']*nSS[i])
            for j in range(nSS[i]):
                
                Density = fetchCenDen_AMReX( j,                  \
                                             FileNumberArray[i], \
                                             DataDirectory[i]    )
                DensityList[i][j] = Density
        SubFrameList = Synchronizer_Density(nDirs,DensityList,nSS)
    
    
    SubFrameLength = len(SubFrameList)
    
    FrameList = [];
    PrevFrame = [0]*nDirs
    FrameList.append(PrevFrame)
    for i in range(SubFrameLength):
        FrameList.append(SubFrameList[i][:])

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
        if gvS.ReferenceBounce:
            MinTime = TimeLists[0][PrevFrame[0]] - BF.BounceTimeList[0]
        else:
            MinTime = TimeLists[0][PrevFrame[0]]
        MinTimeListIndex = [0]
        
#        for i in range(1,nDirs):
        i = 1
        while i < nDirs:
            if gvS.ReferenceBounce:
                CurTime = TimeLists[i][PrevFrame[i]] - BF.BounceTimeList[i]
            else:
                CurTime = TimeLists[i][PrevFrame[i]]
                
#           Locate the next time from all of the lists
            if abs(CurTime - MinTime) <= SyncTol:
#           If more than one list has the minimum time, note it.
                MinTimeListIndex.append(i)
            if CurTime < MinTime-SyncTol:
                MinTime = CurTime
                MinTimeListIndex = [i]
                i = 0

            i += 1
        
#       Advance the frames with minimum time
        for i in range(len(MinTimeListIndex)):
            j = MinTimeListIndex[i]
            PrevFrame[j] = min(PrevFrame[j]+1,nSS[j]-1)
            

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
        SubFrameList = Synchronizer(NewnDirs,NewTimeLists,NewnSS)

 
#       Add returned list to FrameList and fill in frames for finished list(s).
        SubFrameLength = len(SubFrameList)
        for i in range(SubFrameLength):
            for j in range(nDirs):
                if EoL_Flag[j] == 0:
                    PrevFrame[j] = SubFrameList[i][OldToNewMap[j]]+LastFrame[j]
                elif EoL_Flag[j] == 1:
                    PrevFrame[j] = PrevFrame[j]
                
            FrameList.append(PrevFrame[:])

        
#    exit()
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

        MinDen = DensityLists[0][PrevFrame[0]]
        MinDenListIndex = [0]
        
#        for i in range(1,nDirs):
        i = 1
        while i < nDirs:

            CurDen = DensityLists[i][PrevFrame[i]]
                
#           Locate the next time from all of the lists
            if abs(CurDen - MinDen) <= SyncTol:
#           If more than one list has the minimum time, note it.
                MinDenListIndex.append(i)
            if CurDen < MinDen-SyncTol:
                MinDen = CurDen
                MinDenListIndex = [i]
                i = 0

            i += 1
        
#       Advance the frames with minimum time
        for i in range(len(MinDenListIndex)):
            j = MinDenListIndex[i]
            PrevFrame[j] = min(PrevFrame[j]+1,nSS[j]-1)
            

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
                curDen = DensityLists[i][PrevFrame[i]]
                if curDen > gvS.StopDen:
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
        SubFrameList = Synchronizer(NewnDirs,NewDebLists,NewnSS)

 
#       Add returned list to FrameList and fill in frames for finished list(s).
        SubFrameLength = len(SubFrameList)
        for i in range(SubFrameLength):
            for j in range(nDirs):
                if EoL_Flag[j] == 0:
                    PrevFrame[j] = SubFrameList[i][OldToNewMap[j]]+LastFrame[j]
                elif EoL_Flag[j] == 1:
                    PrevFrame[j] = PrevFrame[j]
                
            FrameList.append(PrevFrame[:])

        
#    exit()
    return FrameList
