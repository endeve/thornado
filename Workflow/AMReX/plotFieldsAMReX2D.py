#!/usr/bin/env python3

import yt

# Specify name of problem
ProblemName = 'YahilCollapse_XCFC'

coordinateSystem = 'spherical'

# Specify directory containing amrex Plotfiles
PlotDirectory = '/Users/nickroberts/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/'
#PlotDirectory = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/'

# Specify plot file base name
PlotBaseName = ProblemName + '.plt'

# Specify field to plot
Field = 'PF_V2'

ds = yt.load( PlotDirectory + PlotBaseName + '00000000' )

if ( coordinateSystem == 'spherical' ) :
    axis = 'phi'
else:
    axis = 'z'

yt.SlicePlot( ds, axis, ('boxlib',Field) ).save()
