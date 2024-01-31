#!/usr/bin/env python3



 #=============================================#
#                                               #
#   FindRefinementBoundaries                    #
#                                               #
 #=============================================#
def FindRefinementBoundaries( dX ):

    RefinementLocations = []
    nX = dX.shape[0]
    cur_width = dX[0]
    
    for i in range(nX):
        if ( dX[i] - cur_width != 0.0 ):
            cur_width = dX[i]
            RefinementLocations += [sum(dX[0:i])]


    return RefinementLocations
