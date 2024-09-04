#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
plt.style.use( 'publication.sty' )


import GlobalVariables.Settings     as gvS
import GlobalVariables.Units        as gvU

import Utilities.BounceFinder       as BF
import Utilities.DecadeFinder       as DF
import Utilities.DensityCliffFinder as DCF

from Utilities.GetPlotData          import GetPlotData
from Utilities.Files                import GetFileNumberArray
from Utilities.MakeDataDirectory    import MakeProblemDataDirectory


if __name__ == "__main__":
    #### ========== User Input ==========

    # Specify name of problem
    ProblemName = 'YahilCollapse_XCFC'

    # Specify title of figure
    gvS.FigTitle = ProblemName

    THORNADO_DIR = '/Users/nickroberts/thornado/'
    #THORNADO_DIR = '/home/kkadoogan/Work/Codes/thornado/'

    # Specify directory containing amrex Plotfiles
    gvS.nDirs = 1


    # Specify directory containing amrex Plotfiles
    PlotDirectories = ['None']*gvS.nDirs
    PlotDirectories[0] \
      = THORNADO_DIR + \
          'SandBox/AMReX/Applications/YahilCollapse_XCFC/'

#    PlotDirectories[1] \
#      = THORNADO_DIR + \
#          'SandBox/AMReX/Applications/AdiabaticCollapse_XCFC/AdiabaticCollapse_XCFC_Uni_dr0.5km'
          
          
    gvS.DataType = ['None']*gvS.nDirs
    gvS.DataType[0] = 'AMReX'
#    gvS.DataType[1] = 'AMReX'

    # Specify plot file base name
    PlotBaseName = ProblemName + '.plt'

    # Specify field to plot
    Field = 'PF_D'


    # Specify whether or not to use physical units
    UsePhysicalUnits = True

    # Specify coordinate system (currently supports 'cartesian' and 'spherical')
    CoordinateSystem = 'cartesian'

    # Max level of refinement to plot (-1 plots leaf elements)
    gvS.MaxLevel = -1
    
    # Write extra info to screen
    gvS.Verbose = True

    gvS.ReferenceBounce  = True
    gvS.StopTime         = 9146.03


    #### ====== End of User Input =======
    
 
    DataDirectories = ['None']*gvS.nDirs

    DataDirectories[0] = 'DataDirectories/{:s}_Today'.format( ProblemName )
 
    gvU.SetSpaceTimeUnits(CoordinateSystem, UsePhysicalUnits)
 
#   Append "/" to PlotDirectory, if not present
    for i in range(gvS.nDirs):
        if not PlotDirectories[i][-1] == '/': PlotDirectories[i] += '/'
        if not DataDirectories[i][-1] == '/': DataDirectories[i] += '/'
    
    
        
#   Create Data Directories
    FileNumberArrays = [['None']]*gvS.nDirs
    for i in range(gvS.nDirs):
        FileNumberArrays[i] = GetFileNumberArray( PlotDirectories[i],   \
                                                  PlotBaseName,         \
                                                  -1, -1,  1            )

        MakeProblemDataDirectory( FileNumberArrays[i],\
                                  PlotDirectories[i],  \
                                  PlotBaseName,    \
                                  Field,           \
                                  DataDirectories[i],   \
                                  gvS.DataType[i]   )
 
 
 
 
 
#   Create Bounce Information
    BF.CreateBounceData(  PlotDirectories,  \
                          PlotBaseName,     \
                          DataDirectories,  \
                          gvS.DataType      )
                          
#   Create Density Decade Information
    DF.CreateDecadeData(  PlotDirectories,  \
                          PlotBaseName,     \
                          DataDirectories,  \
                          gvS.DataType      )
                          
                          
    DCF.CreateCliffData( PlotDirectories,   \
                         PlotBaseName,      \
                         DataDirectories,   \
                         FileNumberArrays,  \
                         gvS.DataType       )
 
    import os
    os.system( 'rm -rf __pycache__ ' )


