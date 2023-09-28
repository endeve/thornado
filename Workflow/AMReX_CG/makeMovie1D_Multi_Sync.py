#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )


import GlobalVariables.Settings as gvS
import Utilities.BounceFinder   as BF
from GlobalVariables.Units      import SetSpaceTimeUnits

from Utilities.Files            import GetFileNumberArray
from Utilities.MakeDataArray    import MakeProbelmDataDirectory
from Utilities.MovieMaker_Dual  import MakeMovie
#from Utilities.BounceFinder     import FindBounce, BounceDensityList, BounceFrameList, BounceTimeList


if __name__ == "__main__":

    #### ========== User Input ==========

    # Specify name of problem
    ProblemName = 'AdiabaticCollapse_XCFC'

    # Specify title of figure
    gvS.FigTitle = ProblemName

    # Specify directory containing amrex Plotfiles
    gvS.nDirs = 2
    
    PlotDirectories = ['None']*gvS.nDirs

    PlotDirectories[0] = 'Directory of First Data Set'
    PlotDirectories[1] = 'Directory of Second Data Set'
#    PlotDirectories[1] = 'Directory of Third Data Set. Don't forget to change gvS.nDirs.'

    # Specify plot file base name
    PlotBaseName = ProblemName + '.plt'

    # Specify field to plot
    Field = 'GF_Beta_1'

    # Specify to plot in log-scale
    gvS.UseLogScale_X  = True
    gvS.UseLogScale_Y  = True
    gvS.UseLogScale_2D = False

    # Specify whether or not to use physical units
    UsePhysicalUnits = True

    # Specify coordinate system (currently supports 'cartesian' and 'spherical')
    CoordinateSystem = 'spherical'

    # Only use every <plotEvery> plotfile
    PlotEvery = 1

    # First and last snapshots and number of snapshots to include in movie
    SSi = -1 # -1 -> SSi = 0
    SSf = -1 # -1 -> plotfileArray.shape[0] - 1
    nSS = -1 # -1 -> plotfileArray.shape[0]


    # Max level of refinement to plot (-1 plots leaf elements)
    gvS.MaxLevel = -1

    # Include initial conditions in movie?
    gvS.ShowIC = True

    gvS.PlotMesh = False

    # Write extra info to screen
    gvS.Verbose = True

    # Use custom limts for y-axis (1D) or colorbar (2D)
    gvS.UseCustomLimits = False
    gvS.vmin = 0.70
    gvS.vmax = 1.02

    gvS.MovieRunTime = 10.0 # seconds

    gvS.ReferenceBounce = False
    gvS.StopTime        = 10000

    gvS.ShowRefinement = True
    gvS.RefinementLevels = 9

    gvS.amr = True


    


    #### ====== End of User Input =======

    DataDirectories = ['None']*gvS.nDirs
    DataDirectories[0] = 'DataDirectories/{:s}_StaticMesh'.format( ProblemName )
#    DataDirectories[0] = 'DataDirectories/{:s}_FirstRun'.format( ProblemName )
#    DataDirectories[1] = 'DataDirectories/{:s}_SecondRun'.format( ProblemName )
    DataDirectories[1] = 'DataDirectories/{:s}_ThirdRun'.format( ProblemName )

    ID            = '{:s}_{:s}'.format( ProblemName, Field )
    gvS.MovieName     = 'mov.{:s}.mp4'.format( ID )

    # Append "/" to PlotDirectory, if not present
    for i in range(gvS.nDirs):
        if not PlotDirectories[i][-1] == '/': PlotDirectories[i] += '/'
    
    
    
    #if type(Field) is not list: Field = [ Field ]
    #if type(DataDirectory) is not list: DataDirectory = [ DataDirectory ]




    SetSpaceTimeUnits(CoordinateSystem, UsePhysicalUnits)



                             

    FileNumberArrays = [['None']]*gvS.nDirs
    for i in range(gvS.nDirs):
        FileNumberArrays[i] = GetFileNumberArray( PlotDirectories[i],     \
                                                  PlotBaseName,       \
                                                  SSi, SSf,           \
                                                  PlotEvery           )

        MakeProbelmDataDirectory( FileNumberArrays[i],\
                                  PlotDirectories[i],  \
                                  PlotBaseName,    \
                                  Field,           \
                                  DataDirectories[i]   )



    if gvS.ReferenceBounce:
        
        print( '\n  Finding Bounce Data' )
        print( '  -------------------' )
        BF.BounceDensityList = [0.0]*gvS.nDirs
        BF.BounceTimeList    = [0.0]*gvS.nDirs
        BF.BounceFrameList   = [0]*gvS.nDirs
        for i in range(gvS.nDirs):
            BFrame, BTime, BDensity = BF.FindBounce( PlotDirectories[i],   \
                                                  PlotBaseName          )
            
            BF.BounceDensityList[i] = BDensity
            BF.BounceFrameList[i]   = BFrame
            BF.BounceTimeList[i]    = BTime
    

    MakeMovie( FileNumberArrays, \
               [Field],                              \
               DataDirectories,     \
               Action = 'None'    )





    import os
    os.system( 'rm -rf __pycache__ ' )

