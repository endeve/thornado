#!/usr/bin/env python3

global tb
tb = -1.0
global yLabel
yLabel = ''
global yScale
yScale = 1.0

global ProblemName
global FigTitle
global MovieName

global PlotDirectory
global PlotBaseName

global Field


global UseLogScale_X
global UseLogScale_Y
global UseLogScale_2D

global UsePhysicalUnits

global CoordinateSystem

global PlotEvery

global SSi
global SSf
global nSS


global MaxLevel

global ShowIC
global PlotMesh

global Verbose

global UseCustomLimits
global vmin
global vmax

global MovieRunTime

global ReferenceBounce
global StopTime

global ShowRefinement
global RefinementLevels


global DataDirectory


global nSources

global amr


# Default Settings

ProblemName = ''
FigTitle    = ''

PlotDirectory = ''
PlotBaseName  = ''

Field = ''

UseLogScale_X = True
UseLogScale_Y = True
UseLogScale_2D = True

UsePhysicalUnits = True

CoordinateSystem = 'cartesian'

PlotEvery = 1

SSi = -1
SSf = -1
nSS = -1

MaxLevel = -1

ShowIC   = False
PlotMesh = False

Verbose  = False

UseCustomLimits = False
vmin = 2.0
vmax = 1.0

MovieRunTime = 10

ShowRefinement       = False
RefinementLocations  = []


DataDirectory = ''

nSources = 1

amr = False
