#!/usr/bin/env python3

import numpy as np
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
    ProblemName = 'AdiabaticCollapse_XCFC'

    # Specify title of figure
    gvS.FigTitle = '{:} with AMR'.format(ProblemName)

    THORNADO_DIR = '/Users/nickroberts/thornado/'
    #THORNADO_DIR = '/home/kkadoogan/Work/Codes/thornado/'

    # Specify directory containing amrex Plotfiles
    PlotDirectory \
      = THORNADO_DIR + \
          'SandBox/AdiabaticCollapse_XCFC/Output'
    gvS.DataType = 'Native'

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
    gvS.ShowIC = False

    gvS.PlotMesh = False

    # Write extra info to screen
    gvS.Verbose = True

    # Use custom limts for y-axis
    gvS.UseCustomLimits = False
    gvS.vmin            = -35000.0
    gvS.vmax            = 0.0

    gvS.MovieRunTime = 10.0 # seconds

    gvS.ShowRefinement = True
    gvS.RefinementLevels = 9

    gvS.amr = True

    #### ====== End of User Input =======

    DataDirectory = 'DataDirectories/{:s}'.format( ProblemName )

    ID            = '{:s}_{:s}'.format( ProblemName, Field )
    gvS.MovieName = 'mov.{:s}.mp4'.format( ID )

    # Append "/" to PlotDirectory, if not present
    if not PlotDirectory[-1] == '/': PlotDirectory += '/'

    #if type(Field) is not list: Field = [ Field ]
    #if type(DataDirectory) is not list: DataDirectory = [ DataDirectory ]

    SetSpaceTimeUnits(CoordinateSystem, UsePhysicalUnits)

    FileNumberArray = GetFileNumberArray( PlotDirectory,      \
                                          PlotBaseName,       \
                                          SSi, SSf,           \
                                          PlotEvery           )

    MakeProbelmDataDirectory( FileNumberArray, \
                              PlotDirectory,   \
                              PlotBaseName,    \
                              Field,           \
                              DataDirectory,   \
                              gvS.DataType     )

    MakeMovie( [FileNumberArray], \
               [Field],           \
               [DataDirectory]    )

    import os
    os.system( 'rm -rf __pycache__ ' )
    os.system( 'rm -rf GlobalVariables/__pycache__' )
    os.system( 'rm -rf Utilities/__pycache__' )
