#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )


import GlobalVariables.Settings as gvS
from GlobalVariables.Units   import SetSpaceTimeUnits

from Utilities.Files         import GetFileNumberArray
from Utilities.MakeDataArray import MakeProbelmDataDirectory
from Utilities.MovieMaker    import MakeMovie


if __name__ == "__main__":

    #### ========== User Input ==========

    # Specify name of problem
    ProblemName = 'YahilCollapse_XCFC'

    # Specify title of figure
    gvS.FigTitle = ProblemName

    # Specify directory containing amrex Plotfiles
    PlotDirectoryA = '/Users/nickroberts/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/Data_9Lvls_512/'


    PlotDirectoryB = '/Users/nickroberts/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/Data_9Lvls_512/'

    # Specify plot file base name
    PlotBaseName = ProblemName + '.plt'

    # Specify field to plot
    Field = 'PF_D'

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
    gvS.UseCustomLimits = True
    gvS.vmin = 1.0e-16
    gvS.vmax = 1.0e-1

    gvS.MovieRunTime = 10.0 # seconds

    gvS.ShowRefinement = True
    gvS.RefinementLocations = [ 5.0e+4, 2.5E+4, 1.25E+4, 6.25E+3, 3.125E+3, \
                            1.5625E+3, 7.8125E+2, 3.90625E+2, 1.953125E+2 ]







    #### ====== End of User Input =======

    DataDirectoryA = 'DataDirectories/{:s}_A'.format( ProblemName )
    DataDirectoryB = 'DataDirectories/{:s}_B'.format( ProblemName )

    ID            = '{:s}_{:s}'.format( ProblemName, Field )
    gvS.MovieName     = 'mov.{:s}.mp4'.format( ID )

    # Append "/" to PlotDirectory, if not present
    if not PlotDirectoryA[-1] == '/': PlotDirectoryA += '/'
    if not PlotDirectoryB[-1] == '/': PlotDirectoryB += '/'
    
    
    #if type(Field) is not list: Field = [ Field ]
    #if type(DataDirectory) is not list: DataDirectory = [ DataDirectory ]




    SetSpaceTimeUnits(CoordinateSystem, UsePhysicalUnits)

            
    FileNumberArrayA = GetFileNumberArray( PlotDirectoryA,     \
                                           PlotBaseName,       \
                                           SSi, SSf,           \
                                           PlotEvery           )

    MakeProbelmDataDirectory( FileNumberArrayA,\
                              PlotDirectoryA,  \
                              PlotBaseName,    \
                              Field,           \
                              DataDirectoryA   )


    SSi = 567 # -1 -> SSi = 0
    SSf = 2262
    FileNumberArrayB = GetFileNumberArray( PlotDirectoryB,     \
                                           PlotBaseName,       \
                                           SSi, SSf,           \
                                           PlotEvery           )

    MakeProbelmDataDirectory( FileNumberArrayB,\
                              PlotDirectoryB,  \
                              PlotBaseName,    \
                              Field,           \
                              DataDirectoryB   )



    

    MakeMovie( [FileNumberArrayA, FileNumberArrayB], \
               [Field],                              \
               [DataDirectoryA, DataDirectoryB],     \
               Action = 'RelDiff'    )





    import os
    os.system( 'rm -rf __pycache__ ' )

