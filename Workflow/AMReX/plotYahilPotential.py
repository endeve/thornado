#!/usr/bin/env python3

import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from matplotlib import ticker as mticker
plt.style.use( 'publication.sty' )

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU

from Utilities.YahilProfile import LoadYahilProfile,        \
                                   GetYahilValues,          \
                                   GetYahilPotential,       \
                                   GetYahilMetric,          \
                                   CalcPotential

from Utilities.GetPlotData import GetPlotData
import Utilities.DecadeFinder   as DF

if __name__ == "__main__":

    gamma = 1.3
    central_density = 7e9 # g cm^-3
    pressure = 6e27 # erg cm^-3

    kappa = pressure/central_density**gamma


    YahilFile = '/Users/nickroberts/Files_Needed_For_thornado/YahilHomologousCollapse_Gm_130.dat'

    # Specify name of problem
    ProblemName = 'YahilCollapse_XCFC'

    # Specify title of figure
    FigTitle = ProblemName

    # Specify directory containing amrex Plotfiles
    PlotDirectory = '/Users/nickroberts/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/'

    gvS.DataType = 'AMReX'
    # Specify plot file base name
    PlotBaseName = ProblemName + '.plt'

    # Specify to plot in log-scale
    UseLogScale_X  = True
    UseLogScale_Y  = True
    UseLogScale_2D = False

    # Specify whether or not to use physical units
    UsePhysicalUnits = False

    # Specify coordinate system (currently supports 'cartesian' and 'spherical')
    CoordinateSystem = 'cartesian'

    # Max level of refinement to plot (-1 plots leaf elements)
    MaxLevel = -1

    # Write extra info to screen
    Verbose = True

    # Use custom limts for y-axis (1D) or colorbar (2D)
    UseCustomLimits = False
    vmin = 0.0
    vmax = 2.0

    # Save figure (True) or plot figure (False)
    gvS.SaveFig = False

    # Specify colormap (2D only)
    cmap = 'viridis'

    




    #### ====== End of User Input =======

    polar = False
    if CoordinateSystem == 'spherical':
        polar = True

    # Specify field to plot
    FieldA = 'PF_D'
    FieldB = 'PF_V1'

    FigName = 'fig.YahilTest.png'
    
    colorlist = ['#377eb8', '#ff7f00', '#4daf4a',
                  '#f781bf', '#a65628', '#984ea3',
                  '#999999', '#e41a1c', '#dede00']

    alphalist = [0.5,1.0,1.0]
    widthlist = [3.0,1.0,1.0]
    stylelist = [ '-', '--' ,'-.']
    
    

    

    # Append "/" to PlotDirectory, if not present
    if not PlotDirectory[-1] == '/': PlotDirectory += '/'

    DataDirectory = 'DataDirectories/{:s}_xCFC'.format( ProblemName )

    gvU.SetSpaceTimeUnits(CoordinateSystem, UsePhysicalUnits)



    Decades, DecadeFrames, DecadeTimes = DF.CreateDecadeData( [PlotDirectory],  \
                                                              PlotBaseName,     \
                                                              [DataDirectory],  \
                                                              [gvS.DataType]    )

    print(Decades)
    Frames = DecadeFrames[0]
    Slices = len(Frames)
    print(Frames)
    print(Slices)
    Frames[2] = 2712
    Frames[4] = 4967
#    Frames[5] = 5668
    
    fig, axs = plt.subplots(2,2)
    plt.subplots_adjust(wspace=0, hspace=0)


    
#    Frames[4] = 19867
    for slice in range(Slices):
#    for slice in [1]:
    
        Frame = Frames[slice]
        

        Thor_Density, DensityUnits, Time, \
        X1_C, X2_C, X3_C,                 \
        dX1, dX2, dX3 = GetPlotData( PlotDirectory,     \
                                     PlotBaseName,      \
                                     'PF_D',            \
                                     FrameNumber = Frame,   \
                                     Verbose = Verbose, \
                                     DataType = 'AMReX' )
                                     
        Thor_Velocity, VelocityUnits, Time, \
        X1_C, X2_C, X3_C,                   \
        dX1, dX2, dX3 = GetPlotData( PlotDirectory,     \
                                     PlotBaseName,      \
                                     'PF_V1',           \
                                     FrameNumber = Frame,   \
                                     Verbose = Verbose, \
                                     DataType = 'AMReX' )
                                     
        Thor_Psi, PsiUnits, Time, \
        X1_C, X2_C, X3_C,         \
        dX1, dX2, dX3 = GetPlotData( PlotDirectory,     \
                                     PlotBaseName,      \
                                     'GF_Psi',          \
                                     FrameNumber = Frame,   \
                                     Verbose = Verbose, \
                                     DataType = 'AMReX' )
                                     
        Thor_Alpha, AlphaUnits, Time, \
        X1_C, X2_C, X3_C,             \
        dX1, dX2, dX3 = GetPlotData( PlotDirectory,     \
                                     PlotBaseName,      \
                                     'GF_Alpha',        \
                                     FrameNumber = Frame,   \
                                     Verbose = Verbose, \
                                     DataType = 'AMReX' )

        YahilProfile = LoadYahilProfile( YahilFile )
        
        YahilTimes = [150.0, 51.0, 15.0, 5.0, 1.45, 0.5]
        NewTime = 150-YahilTimes[slice]
        
        Test_Psi,      Test_Alpha     = CalcPotential( X1_C, dX1, Thor_Density)
        Yahil_Density, Yahil_Velocity = GetYahilValues( X1_C, kappa, gamma, NewTime, YahilProfile )
        Yahil_Psi,     Yahil_Alpha    = GetYahilMetric( X1_C, kappa, gamma, NewTime, YahilProfile )
        
        print('Slice',slice)
        print('Frame',Frame)
        print('Time',Time, NewTime)
        print('Cen Den: ',Yahil_Density[0],[Thor_Density[0]])

        
        if (slice == 1 ):
            labellist = ['$\\psi_{h}$','$\\psi_{N}$']
        else:
            labellist = ['_nolegend_','_nolegend_','_nolegend_']
#        labellist = ['$\\rho$','$\\rho_{h}$']

    #   /Density
#        j = 0
#        axs[0,0].plot( X1_C,                        \
#                       Yahil_Density,               \
#                       label     = labellist[j],    \
#                       color     = colorlist[slice],    \
#                       linewidth = widthlist[j],    \
#                       linestyle = stylelist[j],    \
#                       alpha     = alphalist[j]     )
                
        j = 0
        axs[0,0].plot( X1_C,                        \
                       Thor_Density,                \
                       color     = colorlist[slice],    \
                       linewidth = widthlist[j],    \
                       linestyle = stylelist[j],    \
                       alpha     = alphalist[j]     )
         
    #   Velocity
#        j = 0
#        axs[0,1].plot( X1_C,                        \
#                       Yahil_Velocity/299792,       \
#                       color     = colorlist[slice],    \
#                       linewidth = widthlist[j],    \
#                       linestyle = stylelist[j],    \
#                       alpha     = alphalist[j]     )
                
        j = 0
        axs[0,1].plot( X1_C,                        \
                       Thor_Velocity/299792,        \
                       color     = colorlist[slice],    \
                       linewidth = widthlist[j],    \
                       linestyle = stylelist[j],    \
                       alpha     = alphalist[j]     )
         
         
         
    #   Conformal Factor
        j = 0
        axs[1,0].plot( X1_C,                        \
                       Thor_Psi,                    \
                       label     = labellist[j],    \
                       color     = colorlist[slice],    \
                       linewidth = widthlist[j],    \
                       linestyle = stylelist[j],    \
                       alpha     = alphalist[j]     )
                       
        j = 1
        axs[1,0].plot( X1_C,                        \
                       Test_Psi,                   \
                       label     = labellist[j],    \
                       color     = colorlist[slice],    \
                       linewidth = widthlist[j],    \
                       linestyle = stylelist[j],    \
                       alpha     = alphalist[j]     )
                       
#        j = 2
#        axs[1,0].plot( X1_C,                        \
#                       Test_Psi,                    \
#                       label     = labellist[j],    \
#                       color     = colorlist[slice],    \
#                       linewidth = widthlist[j],    \
#                       linestyle = stylelist[j],    \
#                       alpha     = alphalist[j]     )


    #   Lapse Function
        j = 0
        axs[1,1].plot( X1_C,                        \
                       Thor_Alpha,                  \
                       label     = labellist[j],    \
                       color     = colorlist[slice],    \
                       linewidth = widthlist[j],    \
                       linestyle = stylelist[j],    \
                       alpha     = alphalist[j]     )
         
        j = 1
        axs[1,1].plot( X1_C,                        \
                       Test_Alpha,                 \
                       label     = labellist[j],    \
                       color     = colorlist[slice],    \
                       linewidth = widthlist[j],    \
                       linestyle = stylelist[j],    \
                       alpha     = alphalist[j]     )
     
     
     
#   Settings
    xL = 10
    xH = 10**5

    
    axs[0,0].set_xlim( xmin=xL,xmax=xH )
    axs[0,0].set_yscale('log')
    axs[0,0].set_xscale('log')
    axs[0,0].set_ylabel( '$\\rho$' + ' ' + '$'+DensityUnits+'$' )
    
    axs[0,1].set_xlim( xmin=xL,xmax=xH )
    axs[0,1].set_xscale('log')
    axs[0,1].yaxis.set_label_position("right")
    axs[0,1].yaxis.tick_right()
    axs[0,1].set_ylabel( 'v/c' + ' ' + '$'+VelocityUnits+'$' )
    
    axs[1,0].legend()
    axs[1,0].set_xlim( xmin=xL,xmax=xH )
    axs[1,0].set_xscale('log')
    axs[1,0].yaxis.tick_left()
        
    axs[1,1].set_xlim( xmin=xL,xmax=xH )
    axs[1,1].set_xscale('log')
    axs[1,1].yaxis.set_label_position("right")
    axs[1,1].yaxis.tick_right()


    

    xticks = [ 1.0e0, 1.0e1, 1.0e2, 1.0e3, 1.0e4, 1.0e5 ]
    xticklabels = [ r'$10^{0}$', r'$10^{1}$', r'$10^{2}$', r'$10^{3}$',r'$10^{4}$',r'$10^{5}$' ]

    for i in range( axs.shape[0] ):
        for j in range( axs.shape[1] ):
            axs[i,j].set_xlim( 0.0 + 0.25 * 0.5, 8.0e3 + 2.0e3 )
            axs[i,j].set_xscale( 'log' )
            axs[i,j].grid( axis = 'x', which = 'major' )
            axs[i,j].set_xticks( xticks )
            axs[i,j].xaxis.set_minor_locator \
              ( mticker.LogLocator( numticks = 999, subs = 'auto' ) )

    axs[0,0].set_xticklabels( '' )
    axs[0,1].set_xticklabels( '' )
    axs[1,0].set_xticklabels( xticklabels )
    axs[1,1].set_xticklabels( xticklabels )

    if gvS.SaveFig:

        plt.savefig( FigName, dpi = 300 )

    else:

        plt.show()
        
        
    plt.close()


    import os
    os.system( 'rm -rf __pycache__ ' )





