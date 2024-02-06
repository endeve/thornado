#!/usr/bin/env python3

import yt

# Specify name of problem
ProblemName = 'KelvinHelmholtz2D'

# Specify directory containing amrex Plotfiles
#PlotDirectory = '/Users/nickroberts/thornado/SandBox/AMReX/Applications/YahilCollapse_XCFC/'
PlotDirectory = '/home/kkadoogan/Work/Codes/thornado/SandBox/AMReX/dgExperiments_Euler_Relativistic_IDEAL/'

# Specify plot file base name
PlotBaseName = ProblemName + '.plt'

# Specify field to plot
Field = 'PF_D'

ds = yt.load( PlotDirectory + PlotBaseName + '00000000' )

yt.SlicePlot( ds, "z", ("boxlib",Field) ).save()
