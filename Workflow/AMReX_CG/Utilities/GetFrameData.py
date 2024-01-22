#!/usr/bin/env python3

import numpy as np
import h5py as h5
import gc
from charade import detect

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU

from Utilities.Checks       import CoordSystemCheck
from Utilities.CellAverager import computeCellAverage
from Utilities.Mesh         import CreateLocations

#=============================================#
#   Included Routines
#
#   GetFrameData
#   GetFieldData
#
#=============================================#


 #=============================================#
#                                               #
#   GetFrameData                                #
#                                               #
 #=============================================#
def GetFrameData( FilePath,              \
                  DataType,              \
                  Field,                 \
                  SaveTime   = True      ):

    import yt
    import numpy as np

    CoordSystemCheck()

    X1, X2, X3, dX1, dX2, dX3, xL, xH               \
        = CreateLocations( FilePath, TypeIn = DataType )
    
    if DataType.lower() == "amrex":
        # https://yt-project.org/doc/faq/index.html#how-can-i-change-yt-s-log-level
        yt.funcs.mylog.setLevel(40) # Suppress yt warnings

        ds = yt.load( '{:}'.format( FilePath ) )

        Time = ds.current_time.to_ndarray()


        # --- Get Data ---
        Data, DataUnits = GetFieldData( ds,                     \
                                        Field,                  \
                                        gvS.CoordinateSystem,   \
                                        X1, X2, X3,             \
                                        dX1, dX2, dX3           )

        #---- Clean Up ----
        del ds
        gc.collect()
        
    
    elif DataType.lower() == "native":
        Data, DataUnits,    \
        X1, X2, X3,         \
        dX1, dX2, dX3,      \
        Time                \
             = GetFieldData_Native( FilePath,               \
                                    Field,                  \
                                    gvS.CoordinateSystem,   \
                                    X1, X2, X3,             \
                                    dX1, dX2, dX3           )
    
    if not gvS.UsePhysicalUnits: DataUnits = '[]'
    else:                    DataUnits = '[' + DataUnits + ']'


    #---- Clean Up ----
    gc.collect()


    #---- Return ----
    if SaveTime :

        return Data, DataUnits, X1, X2, X3, dX1, dX2, dX3, Time

    else:

        return Data, DataUnits, X1, X2, X3, dX1, dX2, dX3





 #=============================================#
#                                               #
#   GetFieldData                                #
#                                               #
 #=============================================#
def GetFieldData( ds,                   \
                  Field,                \
                  CoordinateSystem,     \
                  X1, X2, X3,           \
                  dX1, dX2, dX3         ):


    nX1 = X1.shape[0]
    nX2 = X2.shape[0]
    nX3 = X3.shape[0]

    Locations = [None]*nX1*nX2*nX3
    for k in range(nX3):
        for j in range(nX2):
            for i in range(nX1):
                Here = k*nX1*nX2    \
                     + j*nX1        \
                     + i
                Locations[Here] = np.array([X1[i],X2[j],X3[k]])





    if Field == 'MPIProcess':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = ''

    elif Field == 'PF_D':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = 'g/cm^3'

    elif Field == 'PF_V1':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = 'km/s'

    elif Field == 'PF_V2':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        if   CoordinateSystem == 'cartesian'  : DataUnits = 'km/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'km/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = '1/s'

    elif Field == 'PF_V3':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        if   CoordinateSystem == 'cartesian'  : DataUnits = 'km/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = '1/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = '1/s'

    elif Field == 'PF_E':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = 'erg/cm^3'

    elif Field == 'PF_Ne':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = '1/cm^3'

    elif Field == 'CF_D':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = 'g/cm^3'

    elif Field == 'CF_S1':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = 'g/cm^2/s'

    elif Field == 'CF_S2':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        if   CoordinateSystem == 'cartesian'  : DataUnits = 'g/cm^2/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'g/cm^2/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'g/cm/s'

    elif Field == 'CF_S3':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        if   CoordinateSystem == 'cartesian'  : DataUnits = 'g/cm^2/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'g/cm/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'g/cm/s'

    elif Field == 'CF_E':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = 'erg/cm^3'

    elif Field == 'CF_Ne':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = '1/cm^3'

    elif Field == 'AF_P':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = 'erg/cm^3'

    elif Field == 'AF_Ye':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = ''

    elif Field == 'AF_T':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = 'K'

    elif Field == 'AF_S':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = 'kb/baryon'

    elif Field == 'AF_Cs':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = 'km/s'

    elif Field == 'GF_Gm_11':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = ''

    elif Field == 'GF_Gm_22':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        if   CoordinateSystem == 'cartesian'  : DataUnits = ''
        elif CoordinateSystem == 'cylindrical': DataUnits = ''
        elif CoordinateSystem == 'spherical'  : DataUnits = 'km^2'

    elif Field == 'GF_Gm_33':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        if   CoordinateSystem == 'cartesian'  : DataUnits = ''
        elif CoordinateSystem == 'cylindrical': DataUnits = 'km^2'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'km^2'

    elif Field == 'GF_K_11':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        if   CoordinateSystem == 'cartesian'  : DataUnits = ''
        elif CoordinateSystem == 'cylindrical': DataUnits = ''
        elif CoordinateSystem == 'spherical'  : DataUnits = ''

    elif Field == 'GF_Psi':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = ''

    elif Field == 'GF_Alpha':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = ''

    elif Field == 'GF_Beta_1':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = 'km/s'

    elif Field == 'DF_TCI':
        Data = np.copy( ds.find_field_values_at_points(("boxlib",Field), Locations ) )
        DataUnits = ''

    # --- Derived Fields ---

    elif Field == 'alphaE':
        alpha = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Alpha'), Locations ) )
        D     = np.copy( ds.find_field_values_at_points(("boxlib",'CF_D'), Locations ) )
        tau   = np.copy( ds.find_field_values_at_points(("boxlib",'CF_E'), Locations ) )

        Data = alpha * ( tau + D )
        DataUnits = 'erg/cm^3'

    elif Field == 'pr4':
        p = np.copy( ds.find_field_values_at_points(("boxlib",'AF_P'), Locations ) )

        Data = np.empty( (nX[0],nX[1],nX[2]), np.float64 )

        for iX1 in range( nX[0] ):
            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):
                    Data[iX1,iX2,iX3] = p[iX1,iX2,iX3] \
                                          * ( X1[iX1,iX2,iX3] * 1.0e5 )**4

        DataUnits = 'erg*cm'

    elif Field == 'RelativisticBernoulliConstant':

        c = 2.99792458e10

        rho   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_D'), Locations ) )
        e     = np.copy( ds.find_field_values_at_points(("boxlib",'PF_E'), Locations ) )
        v1    = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V1'), Locations ) )
        v1    = v1*1.0e5
        v2    = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V2'), Locations ) )
        p     = np.copy( ds.find_field_values_at_points(("boxlib",'AF_P'), Locations ) )
        alpha = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Alpha'), Locations ) )
        Gm11  = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm11'), Locations ) )
        Gm22  = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm22'), Locations ) )
        Gm22  = Gm22*(1.0e5)**2


        VSq = Gm11 * v1**2 + Gm22 * v2**2

        h = c**2 + ( e + p ) / rho
        W = 1.0 / np.sqrt( 1.0 - VSq / c**2 )

        B = alpha * h * W

        Data = B

        DataUnits = 'cm^2/s^2'

    elif Field == 'PolytropicConstant':
        PF_D  = np.copy( ds.find_field_values_at_points(("boxlib",'PF_D'), Locations ) )
        AF_P  = np.copy( ds.find_field_values_at_points(("boxlib",'AF_P'), Locations ) )
        AF_Gm = np.copy( ds.find_field_values_at_points(("boxlib",'AF_Gm'), Locations ) )


        Data  = AF_P / PF_D**AF_Gm

        DataUnits = 'erg/cm^3/(g/cm^3)^(Gamma_{IDEAL})'

    elif Field == 'NonRelativisticSpecificEnthalpy':

        e   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_E'), Locations ) )
        p   = np.copy( ds.find_field_values_at_points(("boxlib",'AF_P'), Locations ) )
        rho = np.copy( ds.find_field_values_at_points(("boxlib",'PF_D'), Locations ) )

        Data = ( e + p ) / rho

        DataUnits = 'cm^2/s^2'

    elif Field == 'RelativisticSpecificEnthalpy':

        c = 2.99792458e10

        e   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_E'), Locations ) )
        p   = np.copy( ds.find_field_values_at_points(("boxlib",'AF_P'), Locations ) )
        rho = np.copy( ds.find_field_values_at_points(("boxlib",'PF_D'), Locations ) )

        Data = ( c**2 + ( e + p ) / rho ) / c**2

        DataUnits = ''

    elif Field == 'LorentzFactor':

        c = 2.99792458e5

        Gm11 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_11'), Locations ) )
        Gm22 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_11'), Locations ) )
        Gm33 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_11'), Locations ) )

        V1   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V1'), Locations ) )
        V2   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V2'), Locations ) )
        V3   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V3'), Locations ) )

        VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

        Data = 1.0 / np.sqrt( 1.0 - VSq / c**2 )

        DataUnits = ''

    elif Field == 'TurbulentVelocity':

        Psi = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Psi'), Locations ) )
        Gm11 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_11'), Locations ) )
        Gm22 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_33'), Locations ) )
        Gm33 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_33'), Locations ) )

        rho = np.copy( ds.find_field_values_at_points(("boxlib",'PF_D'), Locations ) )
        V1   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V1'), Locations ) )
        V2   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V2'), Locations ) )
        V3   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V3'), Locations ) )

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

        Psi = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Psi'), Locations ) )
        Gm11 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_11'), Locations ) )
        Gm22 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_22'), Locations ) )
        Gm33 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_33'), Locations ) )

        rho = np.copy( ds.find_field_values_at_points(("boxlib",'PF_D'), Locations ) )
        V1   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V1'), Locations ) )
        V2   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V2'), Locations ) )
        V3   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V3'), Locations ) )

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


        h1 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_h_1'), Locations ) )
        h2 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_h_2'), Locations ) )
        V1 = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V1'), Locations ) )
        V2 = np.copy( ds.find_field_values_at_points(("boxlib",'AF_V2'), Locations ) )

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

    Data /= gvS.yScale

    return Data, DataUnits









 #=============================================#
#                                               #
#   GetFieldData_Native                         #
#                                               #
 #=============================================#
def GetFieldData_Native( Path,                  \
                         Field,                 \
                         CoordinateSystem,      \
                         X1, X2, X3,            \
                         dX1, dX2, dX3          ):
                  
    gvS.cellAverageFlag = True
         
    PathRoot = Path[:-6]
    FrameNumber = Path[-6:]
    FF_root = PathRoot + 'FluidFields_'
    GF_root = PathRoot + 'GeometryFields_'



    DataFileName_FF = FF_root + str( FrameNumber ) + '.h5'
    DataFileName_GF = GF_root + str( FrameNumber ) + '.h5'
    
    Data_FF      = h5.File( DataFileName_FF, 'r' )
    Data_GF      = h5.File( DataFileName_GF, 'r' )
    
    FF = Data_FF['Fluid Fields']
    GF = Data_GF['Geometry Fields']

    # Second level groups

    PF = FF[ 'Primitive' ]
    CF = FF[ 'Conserved' ]
    AF = FF[ 'Auxiliary' ]
    DF = FF[ 'Diagnostic' ]
    
    
    Time = Data_FF[ 'Time' ][0]
    
    if Field == 'PF_D':
        Data = PF[ 'Comoving Baryon Density'   ][:][0][0]
        if gvS.cellAverageFlag:
            nX1 = len(dX1)
            SqrtGM = GF[ 'Sqrt Spatial Metric Determinant'  ][:][0][0]
            Data, X1 = computeCellAverage(Data, SqrtGM, X1, nX1 )
        DataUnits = 'g/cm^3'

    elif Field == 'PF_V1':
        Data = PF[ 'Three-Velocity (1)'        ][:][0][0]
        if gvS.cellAverageFlag:
            nX1 = len(dX1)
            SqrtGM = GF[ 'Sqrt Spatial Metric Determinant'  ][:][0][0]
            Data, X1 = computeCellAverage(Data, SqrtGM, X1, nX1 )
        DataUnits = 'km/s'

    elif Field == 'PF_V2':
        Data = PF[ 'Three-Velocity (2)'        ][:][0][0]
        if   CoordinateSystem == 'cartesian'  : DataUnits = 'km/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'km/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = '1/s'

    elif Field == 'PF_V3':
        Data = PF[ 'Three-Velocity (3)'        ][:][0][0]
        if   CoordinateSystem == 'cartesian'  : DataUnits = 'km/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = '1/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = '1/s'

    elif Field == 'PF_E':
        Data = PF[ 'Internal Energy Density'   ][:][0][0]
        DataUnits = 'erg/cm^3'

    elif Field == 'PF_Ne':
        Data = PF[ 'Comoving Electron Density' ][:][0][0]
        DataUnits = '1/cm^3'

    elif Field == 'CF_D':
        Data = CF[ 'Conserved Baryon Density'       ][:][0][0]
        DataUnits = 'g/cm^3'

    elif Field == 'CF_S1':
        Data = CF[ 'Conserved Momentum Density (1)' ][:][0][0]
        DataUnits = 'g/cm^2/s'

    elif Field == 'CF_S2':
        Data = CF[ 'Conserved Momentum Density (2)' ][:][0][0]
        if   CoordinateSystem == 'cartesian'  : DataUnits = 'g/cm^2/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'g/cm^2/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'g/cm/s'

    elif Field == 'CF_S3':
        Data = CF[ 'Conserved Momentum Density (3)' ][:][0][0]
        if   CoordinateSystem == 'cartesian'  : DataUnits = 'g/cm^2/s'
        elif CoordinateSystem == 'cylindrical': DataUnits = 'g/cm/s'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'g/cm/s'

    elif Field == 'CF_E':
        Data = CF[ 'Conserved Energy Density'       ][:][0][0]
        DataUnits = 'erg/cm^3'

    elif Field == 'CF_Ne':
        Data = CF[ 'Conserved Electron Density'     ][:][0][0]
        DataUnits = '1/cm^3'

    elif Field == 'AF_P':
        Data = AF[ 'Pressure'                        ][:][0][0]
        DataUnits = 'erg/cm^3'

    elif Field == 'AF_Ye':
        Data = AF[ 'Electron Fraction'               ][:][0][0]
        DataUnits = ''

    elif Field == 'AF_T':
        Data = AF[ 'Temperature'                     ][:][0][0]
        DataUnits = 'K'

    elif Field == 'AF_S':
        Data = AF[ 'Entropy Per Baryon'              ][:][0][0]
        DataUnits = 'kb/baryon'
        if gvS.cellAverageFlag:
            nX1 = len(dX1)
            SqrtGM = GF[ 'Sqrt Spatial Metric Determinant'  ][:][0][0]
            Data, X1 = computeCellAverage(Data, SqrtGM, X1, nX1 )

    elif Field == 'AF_Cs':
        Data = AF[ 'Sound Speed'                     ][:][0][0]
        DataUnits = 'km/s'

    elif Field == 'GF_Gm_11':
        Data = GF[ 'Spatial Metric Component (11)'    ][:][0][0]
        DataUnits = ''

    elif Field == 'GF_Gm_22':
        Data = GF[ 'Spatial Metric Component (22)'    ][:][0][0]
        if   CoordinateSystem == 'cartesian'  : DataUnits = ''
        elif CoordinateSystem == 'cylindrical': DataUnits = ''
        elif CoordinateSystem == 'spherical'  : DataUnits = 'km^2'

    elif Field == 'GF_Gm_33':
        Data = GF[ 'Spatial Metric Component (33)'    ][:][0][0]
        if   CoordinateSystem == 'cartesian'  : DataUnits = ''
        elif CoordinateSystem == 'cylindrical': DataUnits = 'km^2'
        elif CoordinateSystem == 'spherical'  : DataUnits = 'km^2'

    elif Field == 'GF_K_11':
        Data = GF[ 'Extrinsic Curvature Comp. (11)'   ][:][0][0]
        if   CoordinateSystem == 'cartesian'  : DataUnits = ''
        elif CoordinateSystem == 'cylindrical': DataUnits = ''
        elif CoordinateSystem == 'spherical'  : DataUnits = ''

    elif Field == 'GF_Psi':
        Data = GF[ 'Conformal Factor'                 ][:][0][0]
        if gvS.cellAverageFlag:
            nX1 = len(dX1)
            SqrtGM = GF[ 'Sqrt Spatial Metric Determinant'  ][:][0][0]
            Data, X1 = computeCellAverage(Data, SqrtGM, X1, nX1 )
        DataUnits = ''

    elif Field == 'GF_Alpha':
        Data = GF[ 'Lapse Function'                   ][:][0][0]
        if gvS.cellAverageFlag:
            nX1 = len(dX1)
            SqrtGM = GF[ 'Sqrt Spatial Metric Determinant'  ][:][0][0]
            Data, X1 = computeCellAverage(Data, SqrtGM, X1, nX1 )
        DataUnits = ''

    elif Field == 'GF_Beta_1':
        Data = GF[ 'Shift Vector (1)'                 ][:][0][0]
        DataUnits = 'km/s'

    elif Field == 'DF_TCI':
        Data = DF[ 'TCI'        ][:][0][0]
        DataUnits = ''

    # --- Derived Fields ---

    elif Field == 'alphaE':
        alpha = GF[ 'Lapse Function'                   ][:][0][0]
        D     = CF[ 'Conserved Baryon Density'         ][:][0][0]
        tau   = CF[ 'Conserved Energy Density'         ][:][0][0]

        Data = alpha * ( tau + D )
        DataUnits = 'erg/cm^3'

    elif Field == 'pr4':
        p = AF[ 'Pressure'                        ][:][0][0]

        Data = np.empty( (nX[0],nX[1],nX[2]), np.float64 )

        for iX1 in range( nX[0] ):
            for iX2 in range( nX[1] ):
                for iX3 in range( nX[2] ):
                    Data[iX1,iX2,iX3] = p[iX1,iX2,iX3] \
                                          * ( X1[iX1,iX2,iX3] * 1.0e5 )**4

        DataUnits = 'erg*cm'

    elif Field == 'RelativisticBernoulliConstant':

        c = 2.99792458e10
    
        rho   = PF[ 'Comoving Baryon Density'   ][:][0][0]
        e     = PF[ 'Internal Energy Density'   ][:][0][0]
        v1    = PF[ 'Three-Velocity (1)'        ][:][0][0]
        v1    = v1*1.0e5
        v2    = PF[ 'Three-Velocity (2)'        ][:][0][0]
        p     = AF[ 'Pressure'                  ][:][0][0]
        alpha = GF[ 'Lapse Function'                   ][:][0][0]
        Gm11  = GF[ 'Spatial Metric Component (11)'    ][:][0][0]
        Gm22  = GF[ 'Spatial Metric Component (22)'    ][:][0][0]
        Gm22  = Gm22*(1.0e5)**2


        VSq = Gm11 * v1**2 + Gm22 * v2**2

        h = c**2 + ( e + p ) / rho
        W = 1.0 / np.sqrt( 1.0 - VSq / c**2 )

        B = alpha * h * W

        Data = B

        DataUnits = 'cm^2/s^2'

    elif Field == 'PolytropicConstant':
        PF_D  = PF[ 'Comoving Baryon Density'   ][:][0][0]
        AF_P  = AF[ 'Pressure'                  ][:][0][0]
        AF_Gm = AF[ 'Ratio of Specific Heats (Gamma)' ][:][0][0]

        if gvS.cellAverageFlag:
            nX1 = len(dX1)
            SqrtGM = GF[ 'Sqrt Spatial Metric Determinant'  ][:][0][0]
            
#            Data = [[PF_D], [AF_P], [AF_Gm]]
#
#            Data, X1 = computeCellAverage(Data, SqrtGM, X1, nX1 )
            PF_D = computeCellAverage(PF_D, SqrtGM, X1, nX1 )[0]
            AF_P = computeCellAverage(AF_P, SqrtGM, X1, nX1 )[0]
            AF_Gm, X1 = computeCellAverage(AF_Gm, SqrtGM, X1, nX1 )

        Data  = AF_P / PF_D**AF_Gm

        DataUnits = 'erg/cm^3/(g/cm^3)^(Gamma_{IDEAL})'

    elif Field == 'NonRelativisticSpecificEnthalpy':
    
        e   = PF[ 'Internal Energy Density'   ][:][0][0]
        p   = AF[ 'Pressure'                  ][:][0][0]
        rho = PF[ 'Comoving Baryon Density'   ][:][0][0]

        Data = ( e + p ) / rho

        DataUnits = 'cm^2/s^2'

    elif Field == 'RelativisticSpecificEnthalpy':

        c = 2.99792458e10
        
        e   = PF[ 'Internal Energy Density'   ][:][0][0]
        p   = AF[ 'Pressure'                  ][:][0][0]
        rho = PF[ 'Comoving Baryon Density'   ][:][0][0]

        Data = ( c**2 + ( e + p ) / rho ) / c**2

        DataUnits = ''

    elif Field == 'LorentzFactor':

        c = 2.99792458e5

        Gm11 = GF[ 'Spatial Metric Component (11)'    ][:][0][0]
        Gm22 = GF[ 'Spatial Metric Component (22)'    ][:][0][0]
        Gm33 = GF[ 'Spatial Metric Component (33)'    ][:][0][0]
        
        V1   = PF[ 'Three-Velocity (1)'        ][:][0][0]
        V2   = PF[ 'Three-Velocity (2)'        ][:][0][0]
        V3   = PF[ 'Three-Velocity (3)'        ][:][0][0]

        VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

        Data = 1.0 / np.sqrt( 1.0 - VSq / c**2 )

        DataUnits = ''

    elif Field == 'TurbulentVelocity':

        Psi  = GF[ 'Conformal Factor'                 ][:][0][0]
        Gm11 = GF[ 'Spatial Metric Component (11)'    ][:][0][0]
        Gm22 = GF[ 'Spatial Metric Component (22)'    ][:][0][0]
        Gm33 = GF[ 'Spatial Metric Component (33)'    ][:][0][0]

        rho  = PF[ 'Comoving Baryon Density'   ][:][0][0]
        V1   = PF[ 'Three-Velocity (1)'        ][:][0][0]
        V2   = PF[ 'Three-Velocity (2)'        ][:][0][0]
        V3   = PF[ 'Three-Velocity (3)'        ][:][0][0]

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

        Psi  = GF[ 'Conformal Factor'                 ][:][0][0]
        Gm11 = GF[ 'Spatial Metric Component (11)'    ][:][0][0]
        Gm22 = GF[ 'Spatial Metric Component (22)'    ][:][0][0]
        Gm33 = GF[ 'Spatial Metric Component (33)'    ][:][0][0]

        rho  = PF[ 'Comoving Baryon Density'   ][:][0][0]
        V1   = PF[ 'Three-Velocity (1)'        ][:][0][0]
        V2   = PF[ 'Three-Velocity (2)'        ][:][0][0]
        V3   = PF[ 'Three-Velocity (3)'        ][:][0][0]

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


        h1 = GF[ 'Spatial Scale Factor (1)'         ][:][0][0]
        h2 = GF[ 'Spatial Scale Factor (2)'         ][:][0][0]
        V1   = PF[ 'Three-Velocity (1)'        ][:][0][0]
        V2   = PF[ 'Three-Velocity (2)'        ][:][0][0]


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
        
        
    return Data, DataUnits, X1, X2, X3, dX1, dX2, dX3, Time
