#!/usr/bin/env python3

import numpy as np
import gc
from charade import detect

import GlobalVariables.Settings as gvS
import GlobalVariables.Units    as gvU

from Utilities.Checks import CoordSystemCheck

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
                  Field,                 \
                  X1, X2, X3,            \
                  dX1, dX2, dX3,         \
                  SaveTime   = True      ):

    import yt
    import numpy as np

    CoordSystemCheck()

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

    if not gvS.UsePhysicalUnits: DataUnits = '[]'
    else:                    DataUnits = '[' + DataUnits + ']'

    #---- Clean Up ----
    del ds
    gc.collect()
    
    
    #---- Return ----
    if SaveTime :

        return Data, DataUnits, Time

    else:

        return Data, DataUnits





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

        DataUnits = 'erg/cm^3/(g/cm^3)^(Gamma_IDEAL)'

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
        V2   = np.copy( ds.find_field_values_at_points(("boxlib",'AF_V2'), Locations ) )
        V3   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V3'), Locations ) )

        VSq = Gm11 * V1**2 + Gm22 * V2**2 + Gm33 * V3**2

        Data = 1.0 / np.sqrt( 1.0 - VSq / c**2 )

        DataUnits = ''

    elif Field == 'TurbulentVelocity':

        Psi = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Psi'), Locations ) )
        Gm11 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_11'), Locations ) )
        Gm22 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_11'), Locations ) )
        Gm33 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_11'), Locations ) )

        rho = np.copy( ds.find_field_values_at_points(("boxlib",'PF_D'), Locations ) )
        V1   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V1'), Locations ) )
        V2   = np.copy( ds.find_field_values_at_points(("boxlib",'AF_V2'), Locations ) )
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
        Gm22 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_11'), Locations ) )
        Gm33 = np.copy( ds.find_field_values_at_points(("boxlib",'GF_Gm_11'), Locations ) )

        rho = np.copy( ds.find_field_values_at_points(("boxlib",'PF_D'), Locations ) )
        V1   = np.copy( ds.find_field_values_at_points(("boxlib",'PF_V1'), Locations ) )
        V2   = np.copy( ds.find_field_values_at_points(("boxlib",'AF_V2'), Locations ) )
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

    return Data, DataUnits
