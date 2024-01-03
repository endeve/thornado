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
                 DataType,          \
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

    Data, DataUnit,     \
    X1_C, X2_C, X3_C,   \
    dX1, dX2, dX3,      \
    Time = GetFrameData( File,       \
                        DataType,    \
                        Field        )


    return Data, DataUnit, Time, X1_C, X2_C, X3_C, dX1, dX2, dX3
