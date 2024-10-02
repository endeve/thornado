#!/usr/bin/env python3

import numpy as np
import os.path as os
import GlobalVariables.Units    as gvU
import GlobalVariables.Settings as gvS


from Utilities.Files                import Overwrite, CheckForField, CleanField
from Utilities.MakeDataDirectory    import MakeProblemDataDirectory
from Utilities.FetchData            import fetchData_AMReX


 #=============================================#
#                                               #
#   CreateCliffData                             #
#                                               #
 #=============================================#
def CreateCliffData(    PlotDirectories,    \
                        PlotBaseName,       \
                        DataDirectories,    \
                        FileNumberArrays,   \
                        DataTypes,          \
                        forceChoiceF = False, \
                        overwriteF   = True ):

                        
    print( '\n  Creating Cliff Data' )
    print( '------------------------' )



    for i in range(gvS.nDirs):
    
        nFrames = len(FileNumberArrays[i])
    
    
        #   Check to see if any of the existing data directories
        # contain data related to the specified field.
        PathDataDirectory, Flag = CheckForField('DensityCliff',     \
                                                DataDirectories[i], \
                                                FileNumberArrays[i] )


        #   If data related to specified field exists, ask user
        # if they wish it to be overwritten.
        if Flag:
            overwriteF = Overwrite( PathDataDirectory,         \
                                    ForceChoice = forceChoiceF,\
                                    OW = overwriteF            )


        if overwriteF:

            #   If user wishes to overwrite existing data,
            # this cleans the existing data directories.
            CleanField( 'DensityCliff',     \
                        DataDirectories[i],  \
                        FileNumberArrays[i] )




            MakeProblemDataDirectory(   FileNumberArrays[i],\
                                        PlotDirectories[i], \
                                        PlotBaseName,       \
                                        'PF_D',             \
                                        DataDirectories[i], \
                                        DataTypes[i]        )
            
            
            
            for j in range(nFrames):
            
                print('\rWriting Density Cliff Information: {:}/{:}'.format( j+1, nFrames ), end='')
                
                FileNumber = FileNumberArrays[i][j]
                CliffFile = DataDirectories[i] + '/' + '{:}/{:}.dat'.format( FileNumber,'DensityCliff' )

                
                rho, rhoUnits, X1_C, dX1, Time = fetchData_AMReX(j,                     \
                                                                 FileNumberArrays[i],   \
                                                                 DataDirectories[i],    \
                                                                 'PF_D'                 )
                                                                 
                
                DC_Radius, DC_Index = FindDensityCliff( rho, X1_C, dX1 )


                with open( CliffFile, 'w' ) as FileOut:
                    FileOut.write( '# {:}\n'.format( CliffFile  ) )
                    FileOut.write( '# Cliff Radius [{:}]\n'.format( gvU.X1Units ),  )
                    FileOut.write( str(DC_Radius ) + ' ')
                    FileOut.write( str(DC_Index ) + ' ')
                    FileOut.write( '\n' )
                    


                        
 #=============================================#
#                                               #
#   FindDensityCliff                            #
#                                               #
 #=============================================#
def FindDensityCliff( rho, X1_C, dX1 ):


    Radius = 5.0
    
    nLocs = len(X1_C)
    
    drhodx = [0.0]*nLocs
    for i in range(1,nLocs):
        dX = X1_C[i] - X1_C[i-1]
        drho = rho[i] - rho[i-1]
        drhodx[i] = abs(drho/dX)
        
   
   # measure change relative to central density

    Index_Max = max(range(len(drhodx[2:nLocs])), key=drhodx[2:nLocs].__getitem__)
    Radius = X1_C[Index_Max]
    
    
    return Radius, Index_Max






 #=============================================#
#                                               #
#   FetchDensityCliff                           #
#                                               #
 #=============================================#
def FetchDensityCliff( t, FileNumberArray, DataDirectory ):


    if not DataDirectory[-1] == '/': DataDirectory += '/'

    FileDirectory = DataDirectory + str(FileNumberArray[t]) + '/'
    
    DC_File = FileDirectory + '{:}.dat'.format( 'DensityCliff' )
    
    Check = os.isfile(DC_File)
    
    if Check:
        DC_Data = np.loadtxt(DC_File)
        DC_Radius = DC_Data[0]
        DC_Index  = DC_Data[1]
    else:
        print("FetchDensityFile could not find file: ",DC_File)
        exit()
        
    return DC_Radius, DC_Index
