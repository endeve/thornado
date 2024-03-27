#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU
import Utilities.BounceFinder   as BF

from Utilities.FetchData import fetchData_AMReX


 #=============================================#
#                                               #
#   MakeCentralValueOverTimePlot                #
#                                               #
 #=============================================#
def MakeCentralValueOverTimePlot( Field,            \
                                  FileNumberArray,  \
                                  DataDirectories   ):
                            
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

        CreateCentralValueOverTimePlot( Field,              \
                                        FileNumberArray,    \
                                        DataDirectories     )

    else:

        msg =  "\n MakeCentralValueOverTimePlot Error \n"
        msg += "MakeCentralValueOverTimePlot requires the same number of FileNumberArrays \n"
        msg += "as DataDirectories. One FileNumberArray per DataDirectory \n"
        msg += "must be passed into the routine."

        assert(nFiles == nDirs),msg

    return




 #=============================================#
#                                               #
#   CreateCentralValueOverTimePlot              #
#                                               #
 #=============================================#
def CreateCentralValueOverTimePlot( Field,              \
                                    FileNumberArray,    \
                                    DataDirectory       ):


    for i in range(nDirs):
        if DataDirectory[i][-1] != '/': DataDirectory[i] += '/'
        
    
    
    
#    ID      = '{:s}_{:s}'.format( ProblemName, Field )
#    FigName = 'fig.{:s}.png'.format( ID )

    nPlots = 1
    fig,ax = plt.subplots(1, 1)


    
    for i in range(nDirs):
        nSS = len(FileNumberArray[i])
        TimeVals     = ['None']*nSS
        CentralValue = ['None']*nSS
        RoC          = ['None']*nSS
        SRoC         = ['None']*nSS
        
        for j in range(nSS):
            Data, DataUnits, X1_C, dX1, Time = fetchData_AMReX(j, FileNumberArray[i], DataDirectory[i], Field)

            CentralValue[j] = Data[0]
            TimeVals[j] = Time


    ax.set_title(r'Central Value')
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel( Field  + ' ' + r'$\left[\mathrm{{{:}}}\right]$' \
                                  .format( DataUnits[2:-2] ) )

    ax.grid()
    
    ax.plot( TimeVals[1:nSS], CentralValue[1:nSS], 'k' )
    


    ax.set_yscale( 'log' )
#    plt.axhline(y=8192, color='r',linestyle='-')
    plt.axvline(x = BF.BounceTimeList[i],    \
                alpha = 0.5,                \
                linestyle = '--',           \
                label='Bounce',             \
                color = 'red'            )
    
    ax.legend()
    if gvS.SaveFig:

        FigName = 'fig.{:s}'.format( gvS.FigTitle )
        plt.savefig( gvS.FigTitle, dpi = 300 )

    else:

        plt.show()
