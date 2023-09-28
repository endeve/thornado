#!/usr/bin/env python3

import numpy as np

from Utilities.Files        import GetFileNumberArray
from Utilities.GetFrameData import GetFrameData
from Utilities.Mesh         import CreateLocations



global BounceDensityList
global BounceTimeList
global BounceFrameList


 #=============================================#
#                                               #
#   FindBounce                                  #
#                                               #
 #=============================================#
def FindBounce( PlotDirectory,    \
                PlotBaseName ):


    FileNumberArray = GetFileNumberArray( PlotDirectory,    \
                                          PlotBaseName,     \
                                          -1, -1, 1         )

    NumPltFiles = len(FileNumberArray)

    # GetFrameData needs inputs for these values.
    # But, they aren't used when fetching density.
    # So we create theses phonies.
    X1 = np.array([0.0])
    X2 = np.array([0.0])
    X3 = np.array([0.0])
    dX1 = np.array([0.0])
    dX2 = np.array([0.0])
    dX3 = np.array([0.0])
    
    BounceDensity = 0.0
    for i in range(NumPltFiles):

        PlotFileNumber    = FileNumberArray[i]
        PathPlotDirectory = PlotDirectory       \
                          + PlotBaseName        \
                          + '{:}'.format( str(PlotFileNumber).zfill(8) )


#        X1, X2, X3, dX1, dX2, dX3, xL, xH               \
#            = CreateLocations(  PathPlotDirectory,      \
#                                TypeIn = "Leaf"      )


        Density, DensityUnits, Time \
            = GetFrameData( PathPlotDirectory,  \
                            'PF_D',             \
                            X1, X2, X3,         \
                            dX1, dX2, dX3,      \
                            'True'              )


        CurMaxDensity = max(Density)

        if CurMaxDensity > BounceDensity:
            BounceDensity = CurMaxDensity
            BounceFrame   = i
            BounceTime    = Time


    return BounceFrame, BounceTime, BounceDensity

