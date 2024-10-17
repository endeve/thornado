#!/usr/bin/env python3

import numpy as np
from sys import argv

import os
from os.path import isdir, isfile

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU

from Utilities.Files                import Overwrite, CheckForField, CleanField
from Utilities.MakeDataDirectory    import MakeProblemDataDirectory
from Utilities.FetchData            import fetchData_AMReX

 #==============================================#
#                                                #
#   CreateEnclosedMassData                       #
#                                                #
 #==============================================#
def CreateEnclosedMassData( PlotDirectories,    \
                            PlotBaseName,       \
                            DataDirectories,    \
                            DataTypes,           \
                            FileNumberArrays    ):



    print( '\n  Creating Enclosed Mass Data' )
    print( '--------------------------------' )


    for i in range(gvS.nDirs):

        CreateEnclosedMassFiles( PlotDirectories[i],     \
                                 PlotBaseName,           \
                                 DataDirectories[i],     \
                                 DataTypes[i],           \
                                 FileNumberArrays[i]     )









 #==============================================#
#                                                #
#   CreateEnclosedMassFile                       #
#                                                #
 #==============================================#
def CreateEnclosedMassFiles( PlotDirectory,      \
                             PlotBaseName,       \
                             DataDirectory,      \
                             DataType,           \
                             FileNumberArray     ):


    MakeProblemDataDirectory(   FileNumberArray, \
                                PlotDirectory,   \
                                PlotBaseName,    \
                                'PF_D',          \
                                DataDirectory,   \
                                DataType         )
                          

#    MakeProblemDataDirectory(   FileNumberArray, \
#                                PlotDirectory,   \
#                                PlotBaseName,    \
#                                'GF_Psi',        \
#                                DataDirectory,   \
#                                DataType         )


    nFrames = len(FileNumberArray)
    if not DataDirectory[-1] == '/': DataDirectory += '/'
    for i in range(nFrames):
    
        rho, rhoUnits, X1_C, dX1, Time = fetchData_AMReX(i,FileNumberArray,DataDirectory,'PF_D')
#        CF, CFUnits, X1_C, dX1, Time   = fetchData_AMReX(i,FileNumberArray,DataDirectory,'GF_Psi')



        nRLocs = len(dX1)
        ShellMass = [0]*nRLocs
        EnclosedMass =[0]*nRLocs
        for j in range(nRLocs):
        
            R_Inner = X1_C[j] - 0.5*dX1[j]
            R_Outer = X1_C[j] + 0.5*dX1[j]
            
            InnerShellVolume = 4.0/3.0*np.pi*(X1_C[j]**3 - R_Inner**3)
            OuterShellVolume = 4.0/3.0*np.pi*(R_Outer**3 - X1_C[j]**3)
            
            ShellMass[j] = rho[j]*OuterShellVolume
    
            if j == 1:
                EnclosedMass[j] = InnerShellVolume*rho[j]
                
            else:
                EnclosedMass[j] = EnclosedMass[j-1]     \
                                + ShellMass[j-1]        \
                                + InnerShellVolume*rho[j]
                
    
        PlotfileNumber   = str(FileNumberArray[i])
        FileDirectory    = DataDirectory + PlotfileNumber + '/'
        EnclosedMassFile = FileDirectory + '{:}.dat'.format( 'EnclosedMass')
        DataShape = '{:d}'.format(len(EnclosedMass))

        if isdir(FileDirectory):
            os.system( 'rm -rf {:}'.format( EnclosedMassFile ) )

        if gvS.Verbose:
            print( 'Generating data directory: {:} ({:}/{:})'.format \
                ( FileDirectory, i+1, nFrames ) )

    
        if not isfile( EnclosedMassFile ):


            with open( EnclosedMassFile, 'w' ) as FileOut:
    
                FileOut.write( '# {:}\n'              \
                           .format( EnclosedMassFile  ) )
                FileOut.write( '# Array Shape: {:}\n' \
                           .format( DataShape ) )
                FileOut.write( '# Data Units: {:}\n'  \
                           .format( '[g]' ) )
                FileOut.write( '# Min. value: {:.16e}\n' \
                           .format( min(EnclosedMass) ) )
                FileOut.write( '# Max. value: {:.16e}\n' \
                           .format( max(EnclosedMass) ) )

    
                np.savetxt(FileOut, EnclosedMass)






 #==============================================#
#                                                #
#   CreateEnclosedMassFile                       #
#                                                #
 #==============================================#
def CreateEnclosedMassFilesb( PlotDirectory,      \
                             PlotBaseName,       \
                             DataDirectory,      \
                             DataType,           \
                             FileNumberArray     ):



    if DataDirectory[-1] != '/': DataDirectory += '/'
    if PlotDirectory[-1] != '/': PlotDirectory += '/'

    owF = CheckDirectory( DataDirectory,      \
                          FileNumberArray     )

    
    
    
    NumPltFiles = len(FileNumberArray)

    for i in range(NumPltFiles):


        PlotFileNumber = FileNumberArray[i]
        
        if PlotDataType.lower() == 'amrex':
            PathPlotDirectory = PlotDirectory       \
                              + PlotBaseName    \
                              + '{:}'.format( str(PlotFileNumber).zfill(8) )
        else:
            PathPlotDirectory = PlotDirectory       \
                              + PlotBaseName[:-4]   \
                              + '_'     \
                              + '{:}'.format( str(PlotFileNumber).zfill(6) )
                              
                         
        PathDataDirectory = DataDirectory + str(PlotFileNumber) + '/'
        if isdir(PathDataDirectory) and owF:

            PathFieldDirectory = PathDataDirectory         \
                               + '{:}.dat'.format( 'EnclosedMass' )

            os.system( 'rm -rf {:}'.format( PathFieldDirectory ) )



#            if gvS.Verbose:
#                print( 'Generating data directory: {:} ({:}/{:})'.format \
#                    ( PathDataDirectory, i+1, NumPltFiles ) )

        CalcEnclosedMass( PlotDirectory,    \
                          PlotBaseName,     \
                          DataDirectory,    \
                          DataType,         \
                          PlotFileNumber    )
        
        
        
        
        
        
#        WriteEnclosedMass()





# #==============================================#
##                                                #
##   CalcEnclosedMass                             #
##                                                #
# #==============================================#
#def CalcEnclosedMass(   PlotDirectory,    \
#                        PlotBaseName,     \
#                        DataDirectory,    \
#                        DataType,         \
#                        PlotFileNumber    ):
#
                        
        







# #==============================================#
##                                                #
##   CalcEnclosedMass                             #
##                                                #
# #==============================================#
#def WriteEnclosedMass():





 #==============================================#
#                                                #
#   CheckDirectory                               #
#                                                #
 #==============================================#
def CheckDirectory( DataDirectory,              \
                    FileNumberArray,            \
                    forceChoiceD     = False,   \
                    overwriteD       = True,    \
                    forceChoiceF     = False,   \
                    overwriteF       = True     ):


    overwriteD = Overwrite( DataDirectory, ForceChoice = forceChoiceD, OW = overwriteD )
    
    if overwriteD:

        os.system( 'rm -rf {:}'.format( DataDirectory ) )
        os.system(  'mkdir {:}'.format( DataDirectory ) )

    else: # overwriteD == False


        CheckField( 'EnclosedMass',                 \
                    DataDirectory,                  \
                    FileNumberArray,                \
                    forceChoiceF = forceChoiceF,    \
                    overwriteF = overwriteF         )

        CheckField( 'PF_D',                         \
                    DataDirectory,                  \
                    FileNumberArray,                \
                    forceChoiceF = forceChoiceF,    \
                    overwriteF = overwriteF         )
                    
        CheckField( 'GF_CF',                        \
                    DataDirectory,                  \
                    FileNumberArray,                \
                    forceChoiceF = forceChoiceF,    \
                    overwriteF = overwriteF         )


    return overwriteF


 #==============================================#
#                                                #
#   CheckField                                   #
#                                                #
 #==============================================#
def CheckField( Field,                      \
                DataDirectory,              \
                FileNumberArray,            \
                forceChoiceF     = False,   \
                overwriteF       = True     ):



    #   Check to see if any of the existing data directories
    # contain data related to the specified field.
    PathDataDirectory, Flag = CheckForField(Field,          \
                                            DataDirectory,  \
                                            FileNumberArray )


    #   If data related to specified field exists, ask user
    # if they wish it to be overwritten.
    if Flag:
        overwriteF = Overwrite( PathDataDirectory,          \
                                ForceChoice = forceChoiceF, \
                                OW = overwriteF             )


    if overwriteF:

        #   If user wishes to overwrite existing data,
        # this cleans the existing data directories.
        CleanField( Field,              \
                    DataDirectory,      \
                    FileNumberArray     )

