#!/usr/bin/env python3

from Utilities.GetFrameData import GetFrameData
from Utilities.Mesh         import CreateLocations
from Utilities.Files        import GetFileNumberArray, ChoosePlotFile



#=============================================#
#   Included Routines
#
#   GetPlotData
#
#=============================================#


 #=============================================#
#                                               #
#   GetPlotData                                 #
#                                               #
 #=============================================#
def GetPlotData( PlotDirectory,     \
                 PlotBaseName,      \
                 Field,             \
                 FrameNumber = -1,  \
                 argv = ['a'],      \
                 Verbose = "False"  ):

    FileNumberArray = GetFileNumberArray( PlotDirectory, \
                                          PlotBaseName   )
                                       
    if FrameNumber == -1:
        argv = argv
    else:
        argv = ['a',str(FrameNumber)]

    FileName = ChoosePlotFile( FileNumberArray,     \
                               PlotBaseName,        \
                               argv = argv,         \
                               Verbose = Verbose    )


    File = PlotDirectory + FileName

    X1_C, X2_C, X3_C, dX1, dX2, dX3, xL, xH             \
        = CreateLocations(  File, TypeIn = "Leaf"      )

    Data, DataUnit, Time \
        = GetFrameData( File,               \
                        Field,              \
                        X1_C, X2_C, X3_C,   \
                        dX1, dX2, dX3       )


    return Data, DataUnit, Time, X1_C, X2_C, X3_C, dX1, dX2, dX3, xL, xH
