#!/usr/bin/env python3

import numpy as np

def Overwrite( FileOrDirName, ForceChoice = False, OW = False ):

    if ForceChoice: return OW

    from os.path import isfile, isdir

    OW = True

    if ( isfile( FileOrDirName ) or isdir( FileOrDirName ) ):

        YN = input( '{:} exists. overwrite? (Y/N): '.format( FileOrDirName ) )

        if not YN == 'Y':

            print( 'Not overwriting' )
            OW = False

    return OW


def GetFileArray( PlotFileDataDirectory, PlotFileBaseName ):

    from os import listdir

    FileArray \
      = np.sort( np.array( \
          [ file for file in listdir( PlotFileDataDirectory ) ] ) )

    FileList = []

    for iFile in range( FileArray.shape[0] ):

        sFile = FileArray[iFile]

        if( sFile[0:len(PlotFileBaseName)] == PlotFileBaseName \
              and sFile[len(PlotFileBaseName)+1].isdigit() ):
            FileList.append( sFile )

    FileArray = np.array( FileList )

    if not FileArray.shape[0] > 0:

        msg = 'No files found in {:s}\n'.format( PlotFileDataDirectory )
        msg += 'Double check the path\n'

        assert ( FileArray.shape[0] > 0 ), msg

    return FileArray


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
    if ( File[-1] == '/' ): File = File[:-1]

    if Verbose: print( File )

    return File


def GetData( DataDirectory, PlotFileBaseName, Field, \
             CoordinateSystem, UsePhysicalUnits, argv = [ 'a' ], \
             MaxLevel = -1, iX3_CS = 0, \
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

    FileArray = GetFileArray( DataDirectory, PlotFileBaseName )

    File = ChoosePlotFile( FileArray, PlotFileBaseName, argv = argv, \
                           Verbose = Verbose )

    ds = yt.load( '{:}'.format( DataDirectory + File ) )

    if MaxLevel == -1: MaxLevel = ds.index.max_level
    Time = ds.current_time.to_ndarray()
    nX   = ds.domain_dimensions
    xL   = ds.domain_left_edge
    xU   = ds.domain_right_edge

    """
    https://yt-project.org/doc/reference/api/
    yt.data_objects.construction_data_containers.html#yt.data_objects.
    construction_data_containers.YTCoveringGrid
    """
    CoveringGrid \
      = ds.covering_grid \
          ( level           = MaxLevel, \
            left_edge       = xL, \
            dims            = nX * 2**MaxLevel, \
            num_ghost_zones = 0 )

#    ds.force_periodicity()

    nDimsX = 1
    if nX[1] > 1: nDimsX += 1
    if nX[2] > 1: nDimsX += 1

    # --- Get Mesh ---

    xL = xL.to_ndarray()
    xU = xU.to_ndarray()

    X1 = CoveringGrid['X1_C'].to_ndarray()[:,0,0]
    X2 = CoveringGrid['X2_C'].to_ndarray()[0,:,0]
    X3 = CoveringGrid['X3_C'].to_ndarray()[0,0,:]

    dX1 = CoveringGrid['dX1'].to_ndarray()[:,0,0]
    dX2 = CoveringGrid['dX2'].to_ndarray()[0,:,0]
    dX3 = CoveringGrid['dX3'].to_ndarray()[0,0,:]

    if nDimsX < 3:
        X3  = X3 [0] * np.ones( 1, np.float64 )
        dX3 = dX3[0] * np.ones( 1, np.float64 )
    if nDimsX < 2:
        X2  = X2 [0] * np.ones( 1 )
        dX2 = dX2[0] * np.ones( 1 )

    if   Field == 'MPIProcess':

        Data = CoveringGrid['MPIProcess'].to_ndarray()
        DataUnits = ''

    if   Field == 'PF_D':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = 'g/cm**3'

    elif Field == 'PF_V1':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = 'km/s'

    elif Field == 'PF_V2':

        Data = CoveringGrid[Field].to_ndarray()

        if   CoordinateSystem == 'cartesian'  : DataUnits = 'km/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'km/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = '1/s'

    elif Field == 'PF_V3':

        Data = CoveringGrid[Field].to_ndarray()

        if   CoordinateSystem == 'cartesian'  : DataUnits = 'km/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = '1/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = '1/s'

    elif Field == 'PF_E':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = 'erg/cm**3'

    elif Field == 'PF_Ne':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = '1/cm**3'

    elif Field == 'CF_D':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = 'g/cm**3'

    elif Field == 'CF_S1':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = 'g/cm**2/s'

    elif Field == 'CF_S2':

        Data = CoveringGrid[Field].to_ndarray()

        if   CoordinateSystem == 'cartesian'  : DataUnits = 'g/cm**2/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'g/cm**2/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'g/cm/s'

    elif Field == 'CF_S3':

        Data = CoveringGrid[Field].to_ndarray()

        if   CoordinateSystem == 'cartesian'  : DataUnits = 'g/cm**2/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'g/cm/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'g/cm/s'

    elif Field == 'CF_E':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = 'erg/cm**3'

    elif Field == 'CF_Ne':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = '1/cm**3'

    elif Field == 'AF_P':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = 'erg/cm**3'

    elif Field == 'AF_Cs':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = 'km/s'

    elif Field == 'GF_Gm_11':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = ''

    elif Field == 'GF_Gm_22':

        Data = CoveringGrid[Field].to_ndarray()

        if   CoordinateSystem == 'cartesian'  : DataUnits = ''
        elif CoordinateSystem == 'cylindrical': DataUnits = ''
        elif CoordinateSystem == 'spherical'  : DataUnits = 'km**2'

    elif Field == 'GF_Gm_33':

        Data = CoveringGrid[Field].to_ndarray()

        if   CoordinateSystem == 'cartesian'  : DataUnits = ''
        elif CoordinateSystem == 'cylindrical': DataUnits = 'km**2'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'km**2'

    elif Field == 'GF_Psi':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = ''

    elif Field == 'GF_Alpha':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = ''

    elif Field == 'DF_TCI':

        Data = CoveringGrid[Field].to_ndarray()
        DataUnits = ''

    # --- Derived Fields ---

    elif Field == 'pr4':

        p = CoveringGrid['AF_P'].to_ndarray()

        Data = np.empty( (nX[0],nX[1],nX[2]), np.float64 )

        for iX1 in range( nX[0] ):
            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):
                    Data[iX1,iX2,iX3] = p[iX1,iX2,iX3] * ( X1[iX1] * 1.0e5 )**4

        DataUnits = 'erg*cm'

    elif Field == 'RelativisticBernoulliConstant':

        c = 2.99792458e10

        rho   = CoveringGrid['PF_D'    ].to_ndarray()
        e     = CoveringGrid['PF_E'    ].to_ndarray()
        v1    = CoveringGrid['PF_V1'   ].to_ndarray() * 1.0e5
        v2    = CoveringGrid['PF_V2'   ].to_ndarray()
        p     = CoveringGrid['AF_P'    ].to_ndarray()
        alpha = CoveringGrid['GF_Alpha'].to_ndarray()
        Gm11  = CoveringGrid['GF_Gm11' ].to_ndarray()
        Gm22  = CoveringGrid['GF_Gm22' ].to_ndarray() * ( 1.0e5 )**2

        VSq = Gm11 * v1**2 + Gm22 * v2**2

        h = c**2 + ( e + p ) / rho
        W = 1.0 / np.sqrt( 1.0 - VSq / c**2 )

        B = alpha * h * W

        Data = B

        DataUnits = 'cm**2/s**2'

    elif Field == 'PolytropicConstant':

        PF_D  = CoveringGrid['PF_D' ].to_ndarray()
        AF_P  = CoveringGrid['AF_P' ].to_ndarray()
        AF_Gm = CoveringGrid['AF_Gm'].to_ndarray()

        Data  = AF_P / PF_D**AF_Gm

        DataUnits = 'erg/cm**3/(g/cm**3)**(Gamma_IDEAL)'

    elif Field == 'NonRelativisticSpecificEnthalpy':

        e   = CoveringGrid['PF_E'].to_ndarray()
        p   = CoveringGrid['AF_P'].to_ndarray()
        rho = CoveringGrid['PF_D'].to_ndarray()

        Data = ( e + p ) / rho

        DataUnits = 'cm**2/s**2'

    elif Field == 'RelativisticSpecificEnthalpy':

        c = 2.99792458e10

        e   = CoveringGrid['PF_E'].to_ndarray()
        p   = CoveringGrid['AF_P'].to_ndarray()
        rho = CoveringGrid['PF_D'].to_ndarray()

        Data = ( c**2 + ( e + p ) / rho ) / c**2

        DataUnits = ''

    elif Field == 'LorentzFactor':

        c = 2.99792458e5

        Gm11 = CoveringGrid['GF_Gm_11'].to_ndarray()
        Gm22 = CoveringGrid['GF_Gm_22'].to_ndarray()
        Gm33 = CoveringGrid['GF_Gm_33'].to_ndarray()

        V1 = CoveringGrid['PF_V1'].to_ndarray()
        V2 = CoveringGrid['PF_V2'].to_ndarray()
        V3 = CoveringGrid['PF_V3'].to_ndarray()

        VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

        Data = 1.0 / np.sqrt( 1.0 - VSq / c**2 )

        DataUnits = ''

    elif Field == 'TurbulentVelocity':

        Psi  = CoveringGrid['GF_Psi'  ].to_ndarray()
        Gm11 = CoveringGrid['GF_Gm_11'].to_ndarray()
        Gm22 = CoveringGrid['GF_Gm_22'].to_ndarray()
        Gm33 = CoveringGrid['GF_Gm_33'].to_ndarray()

        rho = CoveringGrid['PF_D' ].to_ndarray()
        V1  = CoveringGrid['PF_V1'].to_ndarray()
        V2  = CoveringGrid['PF_V2'].to_ndarray()
        V3  = CoveringGrid['PF_V3'].to_ndarray()

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
                           * np.sin( X2[iX2] ) * dX1[iX1] * dX2[iX2]

                    AngleAveragedRadialVelocity[iX1] \
                      += V1[iX1,iX2,iX3] * rho[iX1,iX2,iX3] \
                           * Psi[iX1,iX2,iX3]**4 \
                           * np.sin( X2[iX2] ) * dX1[iX1] * dX2[iX2]

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

        Psi  = CoveringGrid['GF_Psi'  ].to_ndarray()
        Gm11 = CoveringGrid['GF_Gm_11'].to_ndarray()
        Gm22 = CoveringGrid['GF_Gm_22'].to_ndarray()
        Gm33 = CoveringGrid['GF_Gm_33'].to_ndarray()

        rho = CoveringGrid['PF_D' ].to_ndarray()
        V1  = CoveringGrid['PF_V1'].to_ndarray()
        V2  = CoveringGrid['PF_V2'].to_ndarray()
        V3  = CoveringGrid['PF_V3'].to_ndarray()

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
                           * np.sin( X2[iX2] ) * dX1[iX1] * dX2[iX2]

                    AngleAveragedRadialVelocity[iX1] \
                      += V1[iX1,iX2,iX3] * rho[iX1,iX2,iX3] \
                           * Psi[iX1,iX2,iX3]**4 \
                           * np.sin( X2[iX2] ) * dX1[iX1] * dX2[iX2]

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

        DataUnits = 'erg/cm**3'

    elif Field == 'Vorticity':

        h1 = CoveringGrid['GF_h_1'].to_ndarray()
        h2 = CoveringGrid['GF_h_2'].to_ndarray()
        V1 = CoveringGrid['PF_V1' ].to_ndarray()
        V2 = CoveringGrid['PF_V2' ].to_ndarray()

        h1A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )
        h2A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )
        V1A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )
        V2A = np.empty( (nX[0],nX[1]+2,nX[2]), np.float64 )

        h1A[:,1:-1,:] = h1
        h2A[:,1:-1,:] = h2
        V1A[:,1:-1,:] = V1
        V2A[:,1:-1,:] = V2

        k = 0

        # --- Apply reflecting boundary conditions in theta ---

        for i in range( nX[0] ):

          h1A[i,0,k] = +h1A[i,1,k]
          h2A[i,0,k] = +h2A[i,1,k]
          V1A[i,0,k] = +V1A[i,1,k]
          V2A[i,0,k] = -V2A[i,1,k]

          h1A[i,nX[1]+1,k] = +h1A[i,nX[1],k]
          h2A[i,nX[1]+1,k] = +h2A[i,nX[1],k]
          V1A[i,nX[1]+1,k] = +V1A[i,nX[1],k]
          V2A[i,nX[1]+1,k] = -V2A[i,nX[1],k]

        # --- Compute vorticity in domain using
        #     central differences for derivatives ---

        Data = np.zeros( (nX[0],nX[1],nX[2]), np.float64 )

        k = 0
        for i in range( 1, nX[0] - 1 ):
            for j in range( 1, nX[1] + 1 ):

                Data[i,j-1,k] \
                  = 1.0 / ( h1A[i,j,k]      * h2A[i,j,k] ) \
                    * ( (   h2A[i+1,j,k]**2 * V2A[i+1,j,k] \
                          - h2A[i-1,j,k]**2 * V2A[i-1,j,k] ) \
                          / ( 2.0 * X1[i] ) \
                      - (   h1A[i,j+1,k]**2 * V1A[i,j+1,k] \
                          - h1A[i,j-1,k]**2 * V1A[i,j-1,k] ) \
                          / ( 2.0 * X2[j-1] ) )

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
        print( '  AF_Cs' )
        print( '  GF_Gm_11' )
        print( '  GF_Gm_22' )
        print( '  GF_Gm_33' )
        print( '  GF_Psi' )
        print( '  GF_Alpha' )
        print( '  DF_TCI' )
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

    if nDimsX == 1:

        Data = Data[:,0,0]

    elif nDimsX == 2:

        Data = Data[:,:,0]

    else:

        Data = Data[:,:,iX3_CS]

    if ReturnTime and ReturnMesh:

        return Data, DataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xU, nX, Time

    elif ReturnTime:

        return Data, DataUnits, Time

    elif ReturnMesh:

        return Data, DataUnits, X1, X2, X3, dX1, dX2, dX3, xL, xU, nX

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
