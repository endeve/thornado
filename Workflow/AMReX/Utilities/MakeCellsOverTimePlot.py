#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU

from Utilities.FetchData import fetchData_AMReX, fetchData_Native



 #=============================================#
#                                               #
#   MakeCellsOverTimePlot                       #
#                                               #
 #=============================================#
def MakeCellsOverTimePlot(  FileNumberArray,    \
                            DataDirectories     ):
                            
    global nDirs
    global nLines
    global nFiles


    # Check if FileNumberArray and DataDirectories are lists
    if type(FileNumberArray) is list:
        nFiles = len(FileNumberArray)
    else:
        nFiles = 1


    if type(DataDirectories) is list:
        nDirs = len(DataDirectories)
    else:
        nDirs = 1

        
        
        
        
    if  nFiles == nDirs:

        CreateCellsOverTimePlot(FileNumberArray,    \
                                DataDirectories     )

    else:

        msg =  "\n MakeCellsOverTimePlot Error \n"
        msg += "MakeCellsOverTimePlot requires the same number of FileNumberArrays \n"
        msg += "as DataDirectories. One FileNumberArray per DataDirectory \n"
        msg += "must be passed into the routine."

        assert(nFiles == nDirs),msg

    return




 #=============================================#
#                                               #
#   CreateCellsOverTimePlot                     #
#                                               #
 #=============================================#
def CreateCellsOverTimePlot(FileNumberArray,    \
                            DataDirectory       ):


    for i in range(nDirs):
        if DataDirectory[i][-1] != '/': DataDirectory[i] += '/'
        
    
    
    
#    ID      = '{:s}_{:s}'.format( ProblemName, Field )
#    FigName = 'fig.{:s}.png'.format( ID )

    nPlots = 1
    fig,ax = plt.subplots(1, 1)



    for i in range(nDirs):
        nSS = len(FileNumberArray[i])
        TimeVals = ['None']*nSS
        NumCells = ['None']*nSS
        for j in range(nSS):
            Data, DataUnits,  \
            X1_C, dX1, Time = fetchData_AMReX(j,                 \
                                                 FileNumberArray[i],\
                                                 DataDirectory[i],  \
                                                 'PF_D'              )  # Field Doesn't Matter

            TimeVals[j] = Time
            NumCells[j] = len(X1_C)

        ax.plot( TimeVals, NumCells, 'k' )


    ax.set_title("Radial Elements Over Time")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel("Radial Elements")

    ax.grid()

    if gvS.SaveFig:

        FigName = 'fig.{:s}'.format( gvS.FigTitle )
        plt.savefig( gvS.FigTitle, dpi = 300 )

    else:

        plt.show()
