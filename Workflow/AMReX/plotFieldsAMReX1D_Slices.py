#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )


import GlobalVariables.Settings     as gvS
import GlobalVariables.Units        as gvU

from Utilities.GetPlotData          import GetPlotData
from Utilities.Files                import GetFileNumberArray

"""

Default usage, plots last Plotfile in PlotDirectory:

  $ python3 plotFieldsAMReX.py

Alernate usage, plot specific file in PlotDirectory:

  $ python3 plotFieldsAMReX.py 10

  will plot the *00000010 Plotfile

"""


if __name__ == "__main__":
    #### ========== User Input ==========

    # Specify name of problem
    ProblemName = 'AdiabaticCollapse_XCFC'

    # Specify title of figure
    gvS.FigTitle = ProblemName

    THORNADO_DIR = '/Users/nickroberts/thornado/'
    #THORNADO_DIR = '/home/kkadoogan/Work/Codes/thornado/'

    # Specify directory containing amrex Plotfiles
    PlotDirectory \
      = THORNADO_DIR + \
          'SandBox/AMReX/Applications/AdiabaticCollapse_XCFC/AdiabaticCollapse_XCFC_dr0.25km'
    gvS.DataType = 'AMReX'

    # Specify plot file base name
    PlotBaseName = ProblemName + '.plt'

    # Specify field to plot
    Field = 'PF_D'

    # Specify to plot in log-scale
    gvS.UseLogScale_X  = True
    gvS.UseLogScale_Y  = True
    gvS.UseLogScale_2D = False

    # Specify whether or not to use physical units
    UsePhysicalUnits = False

    # Specify coordinate system (currently supports 'cartesian' and 'spherical')
    CoordinateSystem = 'cartesian'

    # Max level of refinement to plot (-1 plots leaf elements)
    gvS.MaxLevel = -1


    PlotEvery = 2000
    
    # First and last snapshots and number of snapshots to include in movie
    SSi = -1 # -1 -> SSi = 0
    SSf = -1 # -1 -> plotfileArray.shape[0] - 1
    nSS = -1 # -1 -> plotfileArray.shape[0]
    
    # Write extra info to screen
    gvS.Verbose = True

    # Use custom limts for y-axis (1D) or colorbar (2D)
    gvS.UseCustomLimits = False
    gvS.vmin            = 0.0
    gvS.vmax            = 20.0

    gvS.ShowRefinement = False
    gvS.RefinementLevels = 9

    # Save figure (True) or plot figure (False)
    gvS.SaveFig = False

    # Specify colormap (2D only)
    cmap = 'viridis'

    
#### ====== End of User Input =======
    
    DataDirectory = 'DataDirectories/{:s}_AMR025km'.format( ProblemName )

    if not PlotDirectory[-1] == '/': PlotDirectory += '/'


    FileNumberArray = GetFileNumberArray( PlotDirectory,      \
                                          PlotBaseName,       \
                                          SSi, SSf,           \
                                          PlotEvery           )


    #### ====== End of User Input =======

    polar = False
    if CoordinateSystem == 'spherical':
        polar = True

    ID      = '{:s}_{:s}'.format( ProblemName, Field )
    FigName = 'fig.{:s}_Many.png'.format( ID )

    # Append "/" to PlotDirectory, if not present
    if not PlotDirectory[-1] == '/': PlotDirectory += '/'

    gvU.SetSpaceTimeUnits(CoordinateSystem, UsePhysicalUnits)





    fig, ax = plt.subplots( 1, 1 )




    nFiles = FileNumberArray.shape[0]

    for f in range(nFiles):

        Data, DataUnit, Time,   \
        X1_C, X2_C, X3_C,       \
        dX1, dX2, dX3 = GetPlotData( PlotDirectory,     \
                                     PlotBaseName,      \
                                     Field,             \
                                     argv = ['a',str(FileNumberArray[f])],   \
                                     Verbose = gvS.Verbose,     \
                                     DataType = gvS.DataType  )



        nX1 = X1_C.shape[0]
        nX2 = X2_C.shape[0]
        nX3 = X3_C.shape[0]


        nDims = 1
        if nX2 > 1: nDims += 1
        if nX3 > 1: nDims += 1

        if nDims == 1:

            ax.plot( X1_C, Data, 'k.' )
            
        else:

            print('Multi-D not implemented yet.')
            
         
    xL = X1_C[0 ] - 0.5 * dX1[0 ]
    xH = X1_C[-1] + 0.5 * dX1[-1]
         
    if gvS.UseLogScale_X:
        ax.set_xscale( 'log' )
        xL = max( xL, 0.0 + 0.25 * (X1_C[1]-X1_C[0]) )
        
    if gvS.UseLogScale_Y: ax.set_yscale( 'log' )

    if gvS.UseCustomLimits: ax.set_ylim( vmin, vmax )

    if gvS.ShowRefinement:
        bottom, top = plt.ylim()
        ax.plot( (gvS.RefinementLocations[:], gvS.RefinementLocations[:]), \
                 (top, bottom),     \
                 scaley = False,    \
                 color  = 'red',    \
                 zorder = 0,        \
                 alpha  = 0.4       )


    print( xL, xH)

    ax.set_xlim( xL, xH )
    ax.set_xlabel \
      ( r'$x^{{1}}\ \left[\mathrm{{{:}}}\right]$'.format( gvU.X1Units ), \
        fontsize = 15 )
        
    ax.set_ylabel( Field + ' ' + '$'+DataUnit+'$' )
    ax.grid()








    if gvS.SaveFig:

        plt.savefig( FigName, dpi = 300 )

    else:

        plt.show()

    plt.close()

    import os
    os.system( 'rm -rf __pycache__ ' )

