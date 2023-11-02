#!/usr/bin/env python3

import numpy as np
import gc
from charade import detect

def Overwrite( FileOrDirName, ForceChoice = False, OW = False ):

    if ForceChoice: return OW

    from os.path import isfile, isdir

    OW = True

    if ( isfile( FileOrDirName ) or isdir( FileOrDirName ) ):

        if ( isdir( FileOrDirName ) and FileOrDirName[-1] != '/' ):
            FileOrDirName += '/'

        YN = input( '  {:} exists. Overwrite? (y/N): '.format( FileOrDirName ) )

        if YN == 'Y' or YN == 'y' :
            print( '  Overwriting' )
            OW = True
        else:
            print( '  Not overwriting' )
            OW = False

    return OW


def GetFileArray( plotFileDirectory, plotFileBaseName, \
                  SSi = -1, SSf = -1, nSS = -1 ):

    from os import listdir

    fileArray \
      = np.sort( np.array( \
          [ convert( file ) for file in listdir( plotFileDirectory ) ] ) )

    fileList = []
    for iFile in range( fileArray.shape[0] ):

        sFile = fileArray[iFile]

        if( sFile[0:len(plotFileBaseName)] == plotFileBaseName \
              and sFile[len(plotFileBaseName)+1].isdigit() ) :
            fileList.append( sFile )

    fileArray = np.array( fileList )

    if not fileArray.shape[0] > 0:

        msg = '\n>>>No files found.\n'
        msg += '>>>Double check the path: {:}\n'.format( plotFileDirectory )
        msg += '>>>Double check the PlotFileBaseName: {:}\n' \
               .format( plotFileBaseName )
        msg += '>>> Is it plt_ or just plt?\n'

        assert ( fileArray.shape[0] > 0 ), msg

    if SSi < 0: SSi = 0
    if SSf < 0: SSf = fileArray.shape[0] - 1
    if nSS < 0: nSS = fileArray.shape[0]

    plotFileArray = []
    for i in range( nSS ):
        if nSS > 1:
            iSS = SSi + np.int64( ( SSf - SSi ) / ( nSS - 1 ) * i )
        else:
            iSS = 0
        plotFile = str( fileArray[iSS] )
        if plotFile[-1] == '/' :
            plotFileArray.append( plotFile[0:-1] )
        else:
            plotFileArray.append( plotFile )
    plotFileArray = np.array( plotFileArray )

    return plotFileArray
# END GetFileArray

def convert( s ):

    """from https://www.geeksforgeeks.org/python-character-encoding/"""

    # if in the charade instance
    if isinstance(s, str):
        s = s.encode()

    # retrieving the encoding information
    # from the detect() output
    encoding = detect(s)['encoding']

    if encoding == 'utf-8':
        return s.decode()
    else:
        return s.decode(encoding)

def ChoosePlotFile \
      ( FileArray, PlotFileBaseName = 'plt', argv = [ 'a' ], Verbose = False ):

    if len( argv ) == 1:

        # Get last plotfile in directory

        File = FileArray[-1]

    elif len( argv ) == 2:

        if argv[1][0].isalpha():

            File = argv[1]

        else:

            File = PlotFileBaseName + '{:}'.format( argv[1].zfill(8) )

        FileArray = np.array( File )

    else:

        n = len( argv )

        msg = 'len( argv ) must be > 0 and < 3: len( argv ) = {:d}'.format( n )

        arg = ( n > 0 ) & ( n < 3 )
        print( arg )
        assert arg, msg

    # Remove "/" at end of filename, if present
    if File[-1] == '/' : File = np.copy( File[:-1] )

    if Verbose: print( File )

    return File


def GetData( DataDirectory, PlotFileBaseName, Field, \
             CoordinateSystem, UsePhysicalUnits, argv = [ 'a' ], \
             MaxLevel = -1, iX3_CS = 0, \
             SSi = -1, SSf = -1, nSS = -1, \
             ReturnTime = False, ReturnMesh = False, Verbose = False ):

    import yt
    import numpy as np

    msg = 'Invalid choice of CoordinateSystem: {:s}'.format( CoordinateSystem )
    msg += '\n\nValid Choices:\n'
    msg +=   '--------------\n'
    msg +=   'cartesian\n'
    msg +=   'cylindrical\n'
    msg +=   'spherical'

    assert (    CoordinateSystem == 'cartesian' \
             or CoordinateSystem == 'cylindrical' \
             or CoordinateSystem == 'spherical' ), msg

    # https://yt-project.org/doc/faq/index.html#how-can-i-change-yt-s-log-level
    yt.funcs.mylog.setLevel(40) # Suppress yt warnings

    FileArray = GetFileArray( DataDirectory, PlotFileBaseName, \
                              SSi = SSi, SSf = SSf, nSS = nSS )

    File = ChoosePlotFile( FileArray, PlotFileBaseName, argv = argv, \
                           Verbose = Verbose )

    ds = yt.load( '{:}'.format( DataDirectory + File ) )

    if MaxLevel == -1: MaxLevel = ds.index.max_level
    Time = ds.current_time.to_ndarray()
    nX   = ds.domain_dimensions
    xL   = ds.domain_left_edge
    xH   = ds.domain_right_edge

    nDimsX = 1
    if nX[1] > 1: nDimsX += 1
    if nX[2] > 1: nDimsX += 1

    if nDimsX == 1:
        dim_array = [nX[0] * 2**MaxLevel, nX[1], nX[2]]
    elif nDimsX == 2:
        dim_array = [nX[0] * 2**MaxLevel, nX[1] * 2**MaxLevel, nX[2]]
    else:
        dim_array = nX * 2**MaxLevel
        
    
    """
    https://yt-project.org/doc/reference/api/
    yt.data_objects.construction_data_containers.html#yt.data_objects.
    construction_data_containers.YTCoveringGrid
    """
    if nDimsX == 1:
        CoveringGrid \
          = ds.covering_grid \
              ( level           = MaxLevel, \
                left_edge       = xL, \
                dims            = dim_array, \
                num_ghost_zones = 0 )
    else:
        CoveringGrid \
          = ds.covering_grid \
              ( level           = MaxLevel, \
                left_edge       = xL, \
                dims            = dim_array, \
                num_ghost_zones = nX[0] )

    ds.force_periodicity()

    # --- Get Mesh ---
    xL = np.copy( xL.to_ndarray() )
    xH = np.copy( xH.to_ndarray() )

    X1 = np.copy( CoveringGrid['X1_C'].to_ndarray() )
    X2 = np.copy( CoveringGrid['X2_C'].to_ndarray() )
    X3 = np.copy( CoveringGrid['X3_C'].to_ndarray() )

    dX1 = np.copy( CoveringGrid['dX1'].to_ndarray() )
    dX2 = np.copy( CoveringGrid['dX2'].to_ndarray() )
    dX3 = np.copy( CoveringGrid['dX3'].to_ndarray() )


    if   Field == 'MPIProcess':

        Data = np.copy( CoveringGrid['MPIProcess'].to_ndarray() )
        DataUnits = ''

    elif Field == 'PF_D':
        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = 'g/cm^3'
    elif Field == 'PF_V1':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = 'km/s'

    elif Field == 'PF_V2':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )

        if   CoordinateSystem == 'cartesian'  : DataUnits = 'km/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'km/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = '1/s'

    elif Field == 'PF_V3':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )

        if   CoordinateSystem == 'cartesian'  : DataUnits = 'km/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = '1/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = '1/s'

    elif Field == 'PF_E':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = 'erg/cm^3'

    elif Field == 'PF_Ne':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = '1/cm^3'

    elif Field == 'CF_D':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = 'g/cm^3'

    elif Field == 'CF_S1':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = 'g/cm^2/s'

    elif Field == 'CF_S2':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )

        if   CoordinateSystem == 'cartesian'  : DataUnits = 'g/cm^2/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'g/cm^2/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'g/cm/s'

    elif Field == 'CF_S3':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )

        if   CoordinateSystem == 'cartesian'  : DataUnits = 'g/cm^2/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'g/cm/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'g/cm/s'

    elif Field == 'CF_E':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = 'erg/cm^3'

    elif Field == 'CF_Ne':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = '1/cm^3'

    elif Field == 'AF_P':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = 'erg/cm^3'

    elif Field == 'AF_Ye':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = ''

    elif Field == 'AF_T':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = 'K'

    elif Field == 'AF_S':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = 'kb/baryon'

    elif Field == 'AF_Cs':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = 'km/s'

    elif Field == 'GF_Gm_11':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = ''

    elif Field == 'GF_Gm_22':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )

        if   CoordinateSystem == 'cartesian'  : DataUnits = ''
        elif CoordinateSystem == 'cylindrical': DataUnits = ''
        elif CoordinateSystem == 'spherical'  : DataUnits = 'km^2'

    elif Field == 'GF_Gm_33':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )

        if   CoordinateSystem == 'cartesian'  : DataUnits = ''
        elif CoordinateSystem == 'cylindrical': DataUnits = 'km^2'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'km^2'

    elif Field == 'GF_K_11':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )

        if   CoordinateSystem == 'cartesian'  : DataUnits = ''
        elif CoordinateSystem == 'cylindrical': DataUnits = ''
        elif CoordinateSystem == 'spherical'  : DataUnits = ''

    elif Field == 'GF_Psi':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = ''

    elif Field == 'GF_Alpha':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = ''

    elif Field == 'GF_Beta_1':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = 'km/s'

    elif Field == 'DF_TCI':

        Data = np.copy( CoveringGrid[Field].to_ndarray() )
        DataUnits = ''

    # --- Derived Fields ---

    elif Field == 'alphaE':

        alpha = np.copy( CoveringGrid['GF_Alpha'].to_ndarray() )
        D     = np.copy( CoveringGrid['CF_D'].to_ndarray() )
        tau   = np.copy( CoveringGrid['CF_E'].to_ndarray() )

        Data = alpha * ( tau + D )
        DataUnits = 'erg/cm^3'

    elif Field == 'pr4':

        p = np.copy( CoveringGrid['AF_P'].to_ndarray() )

        Data = np.empty( (nX[0],nX[1],nX[2]), np.float64 )

        for iX1 in range( nX[0] ):
            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):
                    Data[iX1,iX2,iX3] = p[iX1,iX2,iX3] \
                                          * ( X1[iX1,iX2,iX3] * 1.0e5 )**4

        DataUnits = 'erg*cm'

    elif Field == 'RelativisticBernoulliConstant':

        c = 2.99792458e10

        rho   = np.copy( CoveringGrid['PF_D'    ].to_ndarray() )
        e     = np.copy( CoveringGrid['PF_E'    ].to_ndarray() )
        v1    = np.copy( CoveringGrid['PF_V1'   ].to_ndarray() ) * 1.0e5
        v2    = np.copy( CoveringGrid['PF_V2'   ].to_ndarray() )
        p     = np.copy( CoveringGrid['AF_P'    ].to_ndarray() )
        alpha = np.copy( CoveringGrid['GF_Alpha'].to_ndarray() )
        Gm11  = np.copy( CoveringGrid['GF_Gm11' ].to_ndarray() )
        Gm22  = np.copy( CoveringGrid['GF_Gm22' ].to_ndarray() ) * ( 1.0e5 )**2

        VSq = Gm11 * v1**2 + Gm22 * v2**2

        h = c**2 + ( e + p ) / rho
        W = 1.0 / np.sqrt( 1.0 - VSq / c**2 )

        B = alpha * h * W

        Data = B

        DataUnits = 'cm^2/s^2'

    elif Field == 'PolytropicConstant':

        PF_D  = np.copy( CoveringGrid['PF_D' ].to_ndarray() )
        AF_P  = np.copy( CoveringGrid['AF_P' ].to_ndarray() )
        AF_Gm = np.copy( CoveringGrid['AF_Gm'].to_ndarray() )

        Data  = AF_P / PF_D**AF_Gm

        DataUnits = 'erg/cm^3/(g/cm^3)^(Gamma_IDEAL)'

    elif Field == 'NonRelativisticSpecificEnthalpy':

        e   = np.copy( CoveringGrid['PF_E'].to_ndarray() )
        p   = np.copy( CoveringGrid['AF_P'].to_ndarray() )
        rho = np.copy( CoveringGrid['PF_D'].to_ndarray() )

        Data = ( e + p ) / rho

        DataUnits = 'cm^2/s^2'

    elif Field == 'RelativisticSpecificEnthalpy':

        c = 2.99792458e10

        e   = np.copy( CoveringGrid['PF_E'].to_ndarray() )
        p   = np.copy( CoveringGrid['AF_P'].to_ndarray() )
        rho = np.copy( CoveringGrid['PF_D'].to_ndarray() )

        Data = ( c**2 + ( e + p ) / rho ) / c**2

        DataUnits = ''

    elif Field == 'LorentzFactor':

        c = 2.99792458e5

        Gm11 = np.copy( CoveringGrid['GF_Gm_11'].to_ndarray() )
        Gm22 = np.copy( CoveringGrid['GF_Gm_22'].to_ndarray() )
        Gm33 = np.copy( CoveringGrid['GF_Gm_33'].to_ndarray() )

        V1 = np.copy( CoveringGrid['PF_V1'].to_ndarray() )
        V2 = np.copy( CoveringGrid['PF_V2'].to_ndarray() )
        V3 = np.copy( CoveringGrid['PF_V3'].to_ndarray() )

        VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

        Data = 1.0 / np.sqrt( 1.0 - VSq / c**2 )

        DataUnits = ''

    elif Field == 'TurbulentVelocity':

        Psi  = np.copy( CoveringGrid['GF_Psi'  ].to_ndarray() )
        Gm11 = np.copy( CoveringGrid['GF_Gm_11'].to_ndarray() )
        Gm22 = np.copy( CoveringGrid['GF_Gm_22'].to_ndarray() )
        Gm33 = np.copy( CoveringGrid['GF_Gm_33'].to_ndarray() )

        rho = np.copy( CoveringGrid['PF_D' ].to_ndarray() )
        V1  = np.copy( CoveringGrid['PF_V1'].to_ndarray() )
        V2  = np.copy( CoveringGrid['PF_V2'].to_ndarray() )
        V3  = np.copy( CoveringGrid['PF_V3'].to_ndarray() )

        # --- Compute angle-averaged and
        #     mass density weighted radial velocity ---

        AngleAveragedMass           = np.zeros( (nX[0]), np.float64 )
        AngleAveragedRadialVelocity = np.zeros( (nX[0]), np.float64 )

        Data = np.empty( nX, np.float64 )

        for iX1 in range( nX[0] ):

            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):

                    AngleAveragedMass[iX1] \
                      += rho[iX1,iX2,iX3] \
                           * Psi[iX1,iX2,iX3]**4 \
                           * np.sin( X2[iX1,iX2,iX3] ) \
                           * dX1[iX1,iX2,iX3] * dX2[iX1,iX2,iX3]

                    AngleAveragedRadialVelocity[iX1] \
                      += V1[iX1,iX2,iX3] * rho[iX1,iX2,iX3] \
                           * Psi[iX1,iX2,iX3]**4 \
                           * np.sin( X2[iX1,iX2,iX3] ) \
                           * dX1[iX1,iX2,iX3] * dX2[iX1,iX2,iX3]

            AngleAveragedRadialVelocity[iX1] /= AngleAveragedMass[iX1]

            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):

                    Data[iX1,iX2,iX3] \
                      = np.sqrt( \
                          Gm11[iX1,iX2,iX3] \
                            * ( V1[iX1,iX2,iX3] \
                                  - AngleAveragedRadialVelocity[iX1] )**2 \
                            + Gm22[iX1,iX2,iX3] * V2[iX1,iX2,iX3]**2 \
                            + Gm33[iX1,iX2,iX3] * V3[iX1,iX2,iX3]**2 )

        DataUnits = 'km/s'

    elif Field == 'TurbulentEnergyDensity':

        Psi  = np.copy( CoveringGrid['GF_Psi'  ].to_ndarray() )
        Gm11 = np.copy( CoveringGrid['GF_Gm_11'].to_ndarray() )
        Gm22 = np.copy( CoveringGrid['GF_Gm_22'].to_ndarray() )
        Gm33 = np.copy( CoveringGrid['GF_Gm_33'].to_ndarray() )

        rho = np.copy( CoveringGrid['PF_D' ].to_ndarray() )
        V1  = np.copy( CoveringGrid['PF_V1'].to_ndarray() )
        V2  = np.copy( CoveringGrid['PF_V2'].to_ndarray() )
        V3  = np.copy( CoveringGrid['PF_V3'].to_ndarray() )

        AngleAveragedMass           = np.zeros( (nX[0]), np.float64 )
        AngleAveragedRadialVelocity = np.zeros( (nX[0]), np.float64 )

        c = 2.99792458e5

        Data = np.empty( nX, np.float64 )

        for iX1 in range( nX[0] ):

            # --- Compute angle-averaged and
            #     mass density weighted radial velocity ---

            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):

                    AngleAveragedMass[iX1] \
                      += rho[iX1,iX2,iX3] \
                           * Psi[iX1,iX2,iX3]**4 \
                           * np.sin( X2[iX1,iX2,iX3] ) \
                           * dX1[iX1,iX2,iX3] * dX2[iX1,iX2,iX3]

                    AngleAveragedRadialVelocity[iX1] \
                      += V1[iX1,iX2,iX3] * rho[iX1,iX2,iX3] \
                           * Psi[iX1,iX2,iX3]**4 \
                           * np.sin( X2[iX1,iX2,iX3] ) \
                           * dX1[iX1,iX2,iX3] * dX2[iX1,iX2,iX3]

            AngleAveragedRadialVelocity[iX1] /= AngleAveragedMass[iX1]

            # --- Compute turbulent energy density ---

            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):

                    # --- BetaSq = v_i * v^i / c^2 ---

                    BetaSq = ( Gm11[iX1,iX2,iX3] \
                                 * ( V1[iX1,iX2,iX3] \
                                       - AngleAveragedRadialVelocity[iX1] )**2 \
                                 + Gm22[iX1,iX2,iX3] * V2[iX1,iX2,iX3]**2 \
                                 + Gm33[iX1,iX2,iX3] * V3[iX1,iX2,iX3]**2 ) \
                               / c**2

                    W = 1.0 / np.sqrt( 1.0 - BetaSq )

                    Data[iX1,iX2,iX3] \
                      = rho[iX1,iX2,iX3] * ( c * 1.0e5 )**2 \
                          * W**2 * BetaSq / ( W + 1.0 )

        DataUnits = 'erg/cm^3'

    elif Field == 'Vorticity':

        h1 = np.copy( CoveringGrid['GF_h_1'].to_ndarray() )
        h2 = np.copy( CoveringGrid['GF_h_2'].to_ndarray() )
        V1 = np.copy( CoveringGrid['PF_V1' ].to_ndarray() )
        V2 = np.copy( CoveringGrid['PF_V2' ].to_ndarray() )

        h1A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )
        h2A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )
        V1A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )
        V2A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )

        h1A[:,1:-1,:] = np.copy( h1 )
        h2A[:,1:-1,:] = np.copy( h2 )
        V1A[:,1:-1,:] = np.copy( V1 )
        V2A[:,1:-1,:] = np.copy( V2 )

        # --- Apply reflecting boundary conditions in theta ---

        for i in range( nX[0] ):
            for k in range( nX[2] ):

                h1A[i,0,k] = +h1A[i,1,k]
                h2A[i,0,k] = +h2A[i,1,k]
                V1A[i,0,k] = +V1A[i,1,k]
                V2A[i,0,k] = -V2A[i,1,k]

                h1A[i,-1,k] = +h1A[i,-2,k]
                h2A[i,-1,k] = +h2A[i,-2,k]
                V1A[i,-1,k] = +V1A[i,-2,k]
                V2A[i,-1,k] = -V2A[i,-2,k]

        # --- Compute vorticity in domain using
        #     central differences for derivatives (assume 2D) ---

        Data = np.zeros( (nX[0],nX[1],nX[2]), np.float64 )

        for i in range( 1, nX[0] - 1 ):
            for j in range( 1, nX[1] + 1 ):
                for k in range( nX[2] ):

                    Data[i,j-1,k] \
                      = 1.0 / ( h1A[i,j,k]      * h2A[i,j,k] ) \
                        * ( (   h2A[i+1,j,k]**2 * V2A[i+1,j,k] \
                              - h2A[i-1,j,k]**2 * V2A[i-1,j,k] ) \
                              / ( 2.0 * X1[i,j,k] ) \
                          - (   h1A[i,j+1,k]**2 * V1A[i,j+1,k] \
                              - h1A[i,j-1,k]**2 * V1A[i,j-1,k] ) \
                              / ( 2.0 * X2[i,j-1,k] ) )

        DataUnits = '1/s'

    else:

        print( '\nInvalid field: {:}'.format( Field ) )
        print( '\nValid choices:' )
        print( '--------------' )
        print( '  MPIProcess' )
        print( '  PF_D' )
        print( '  PF_V1' )
        print( '  PF_V2' )
        print( '  PF_V3' )
        print( '  PF_E' )
        print( '  PF_Ne' )
        print( '  CF_D' )
        print( '  CF_S1' )
        print( '  CF_S2' )
        print( '  CF_S3' )
        print( '  CF_E' )
        print( '  CF_Ne' )
        print( '  AF_P' )
        print( '  AF_Ye' )
        print( '  AF_T' )
        print( '  AF_S' )
        print( '  AF_Cs' )
        print( '  GF_Gm_11' )
        print( '  GF_Gm_22' )
        print( '  GF_Gm_33' )
        print( '  GF_K_11' )
        print( '  GF_Psi' )
        print( '  GF_Alpha' )
        print( '  GF_Beta_1' )
        print( '  DF_TCI' )
        print( '  alphaE' )
        print( '  pr4' )
        print( '  RelativisticBernoulliConstant' )
        print( '  PolytropicConstant' )
        print( '  NonRelativisticSpecificEnthalpy' )
        print( '  RelativisticSpecificEnthalpy' )
        print( '  LorentzFactor' )
        print( '  TurbulentVelocity' )
        print( '  TurbulentEnergyDensity' )
        print( '  Vorticity' )

        assert 0, 'Invalid choice of field'



    if not UsePhysicalUnits: DataUnits = '[]'
    else:                    DataUnits = '[' + DataUnits + ']'

    del ds, CoveringGrid
    gc.collect()

    if ReturnTime and ReturnMesh:

        return Data, DataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xH, nX, Time

    elif ReturnTime:

        return Data, DataUnits, Time

    elif ReturnMesh:

        return Data, DataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xH, nX

    else:

        return Data, DataUnits


def GetNorm( UseLogScale, Data, vmin = +1.0e100, vmax = -1.0e100 ):

    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, SymLogNorm

    if vmin > +1.0e99: vmin = Data.min()
    if vmax < -1.0e99: vmax = Data.max()

    if UseLogScale:

        if np.any( Data <= 0.0 ):

            Norm = SymLogNorm( vmin = vmin, vmax = vmax, \
                               linthresh = 1.0e-2, base = 10 )

        else:

            Norm = LogNorm   ( vmin = vmin, vmax = vmax )

    else:

        Norm = plt.Normalize ( vmin = vmin, vmax = vmax )

    return Norm

def MapCenterToCorners( X1_C, X2_C, dX1, dX2 ):

    if( X1_C.ndim == 3 ): X1_C = np.copy( X1_C[:,:,0] )
    if( X2_C.ndim == 3 ): X2_C = np.copy( X2_C[:,:,0] )
    if( dX1 .ndim == 3 ): dX1  = np.copy( dX1 [:,:,0] )
    if( dX2 .ndim == 3 ): dX2  = np.copy( dX2 [:,:,0] )

    nX = [ X1_C.shape[0], X1_C.shape[1] ]

    X1c = np.empty( (nX[0]+1,nX[1]+1), np.float64 )
    X2c = np.empty( (nX[0]+1,nX[1]+1), np.float64 )
    for iX1 in range( nX[0] ):
        for iX2 in range( nX[1] ):
            X1c[iX1,iX2] = X1_C[iX1,iX2] - 0.5 * dX1[iX1,iX2]
            X2c[iX1,iX2] = X2_C[iX1,iX2] - 0.5 * dX2[iX1,iX2]
            if   iX2 == nX[1]-1 and iX1 == nX[0]-1:
                X1c[iX1,iX2+1  ] = X1_C[iX1,iX2] - 0.5 * dX1[iX1,iX2]
                X2c[iX1,iX2+1  ] = X2_C[iX1,iX2] + 0.5 * dX2[iX1,iX2]
                X1c[iX1+1,iX2  ] = X1_C[iX1,iX2] + 0.5 * dX1[iX1,iX2]
                X2c[iX1+1,iX2  ] = X2_C[iX1,iX2] - 0.5 * dX2[iX1,iX2]
                X1c[iX1+1,iX2+1] = X1_C[iX1,iX2] + 0.5 * dX1[iX1,iX2]
                X2c[iX1+1,iX2+1] = X2_C[iX1,iX2] + 0.5 * dX2[iX1,iX2]
            elif iX2 == nX[1]-1:
                X1c[iX1,iX2+1] = X1_C[iX1,iX2] - 0.5 * dX1[iX1,iX2]
                X2c[iX1,iX2+1] = X2_C[iX1,iX2] + 0.5 * dX2[iX1,iX2]
            elif iX1 == nX[0]-1:
                X1c[iX1+1,iX2] = X1_C[iX1,iX2] + 0.5 * dX1[iX1,iX2]
                X2c[iX1+1,iX2] = X2_C[iX1,iX2] - 0.5 * dX2[iX1,iX2]

    return X1c, X2c
