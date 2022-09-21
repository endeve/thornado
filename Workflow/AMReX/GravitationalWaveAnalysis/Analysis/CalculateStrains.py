# --- Calculating GW Strains ---

from math import pi

import gc
import numpy as np
import matplotlib.pyplot as plt
import subprocess
import yt
import math as m
from UtilitiesModule import ChoosePlotFile

from Meridional import meridional_int_module
from Meridional import basisf

half  = 1.0 / 2.0
twopi = 2.0 * pi
pio2  = half * pi

# --- Get data from a single file ---

def GetDataR( DataDirectory, PlotFileBaseName, File, FileArray, \
              CoordinateSystem, UsePhysicalUnits, argv = [ 'a' ], \
              Verbose = True ):

    import yt
    import numpy as np

    if Verbose:

       print( '\n  Calling GetData...' )
       print(   '  ------------------' )
       print(   '    DataDirectory:    {:s}'.format( DataDirectory ) )
       print(   '    PlotFileBaseName: {:s}'.format( PlotFileBaseName ) )
       print(   '    Field:            {:s}'.format( Field ) )
       print(   '    CoordinateSystem: {:s}'.format( CoordinateSystem ) )

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

    if Verbose:
       print(   '    File:             {:}\n'.format( File ) )

    ds = yt.load( '{:}'.format( DataDirectory + File ) )

    MaxLevel = ds.index.max_level
    Time     = ds.current_time.to_ndarray()
    nX       = ds.domain_dimensions
    xL       = ds.domain_left_edge
    xU       = ds.domain_right_edge

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
            num_ghost_zones = nX[0] )

    ds.force_periodicity()

    nDimsX = 1
    if nX[1] > 1: nDimsX += 1
    if nX[2] > 1: nDimsX += 1

    # --- Get Mesh ---

    KmToCm = 1.0e5

    xL = xL.to_ndarray()
    xU = xU.to_ndarray()

    dX = ( xU - xL ) / np.float64( nX )

    dX1 = dX[0]
    dX2 = dX[1]
    dX3 = dX[2]

    X1_C = np.linspace( xL[0] + dX[0] / 2.0, xU[0] - dX[0] / 2.0, nX[0] ) * KmToCm
    X2_C = np.linspace( xL[1] + dX[1] / 2.0, xU[1] - dX[1] / 2.0, nX[1] )
    X3_C = np.linspace( xL[2] + dX[2] / 2.0, xU[2] - dX[2] / 2.0, nX[2] )

    X1_F = np.linspace( xL[0], xU[0], nX[0] + 1 ) * KmToCm
    X2_F = np.linspace( xL[1], xU[1], nX[1] + 1 )
    X3_F = np.linspace( xL[2], xU[2], nX[2] + 1 )

    #if File == FileArray[0]:
       #print('X1_C')
       #print(X1_C)
       #print('X2_C')
       #print(X2_C)
       #print('X3_C')
       #print(X3_C)
       #print('X1_F')
       #print(X1_F)
       #print('X2_F')
       #print(X2_F)
       #print('X3_F')
       #print(X3_F)

    D = CoveringGrid['PF_D'].to_ndarray()
    DataUnit_D  = 'g/cm**3'
    D = D[:,:,0]

   #print('Density')
   #print(D)

    V1_L = CoveringGrid['PF_V1'].to_ndarray()
    V1_L = V1_L * KmToCm
    DataUnit_V1 = 'cm/s'
    V1_L = V1_L[:,:,0]

   #print('Velocity 1')
   #print(V1_L)

    V2_A = CoveringGrid['PF_V2'].to_ndarray()
    H2 = CoveringGrid['GF_h_2'].to_ndarray()
    V2_L = H2 * V2_A * KmToCm
    DataUnit_V2 = 'cm/s'
    V2_L = V2_L[:,:,0]

   #print('Velocity 2')
   #print(V2_L)

    V3_A = CoveringGrid['PF_V3'].to_ndarray()
    H3 = CoveringGrid['GF_h_3'].to_ndarray()
    V3_L = H3 * V3_A * KmToCm
    DataUnit_V3 = 'cm/s'
    V3_L = V3_L[:,:,0]

   #print('Velocity 3')
   #print(V3_L)

    if not UsePhysicalUnits: DataUnit = ''

    ds.close()
    del ds, CoveringGrid
    gc.collect()
    return D, V1_L, V2_L, V3_L, X1_C, X2_C, X3_C, \
           X1_F, X2_F, X3_F, xL, xU, nX, dX, Time

# --- Get data Analytically for a single time step ---

def GetData_Analytic_Strain( cvel, Gc, pi, M0, ms, R0, RT0, W0, D0, amp, nX, t ):
    
    X1_F = np.zeros( nX[0]+1 )
    X2_F = np.zeros( nX[1]+1 )
    X3_F = np.zeros( nX[2]+1 )
    X1_C = np.zeros( nX[0] )
    X2_C = np.zeros( nX[1] )
    X3_C = np.zeros( nX[2] )
    dX = np.zeros( 3 )

    D = np.zeros( (nX[0], nX[1]) )
    V1_L = np.zeros( (nX[0], nX[1]) )
    V2_L = np.zeros( (nX[0], nX[1]) )
    V3_L = np.zeros( (nX[0], nX[1]) )

    R = R0 * ( 1.0 + amp*m.cos(W0*t) )
    RT = m.sqrt( ms / (2.0*pi**2*D0*R) )
    
    # Getting radial coordinate geometry

    dX[0] = ( R + RT ) * 1.2 / nX[0]
    
    #for i in range(nX[0]):
    #    X1_F[i] = (R + RT) * 1.2 * i / nX[0]
    #    if i == 0:
    #        X1_C[0] = 0.5 * dX[0]
    #    else:
    #        X1_C[i] = X1_C[i-1] + dX[0]

    #X1_F[nX[0]] = ( R + RT ) * 1.2 

    # Getting theta coordinate geometry
    dX[1] = pi / nX[1]
    dX[2] = 1.0
    
    X1_C = np.linspace( 0.0 + dX[0] / 2.0, 1.2 * ( R + RT ) - dX[0] / 2.0, nX[0] )
    X2_C = np.linspace( 0.0 + dX[1] / 2.0,  pi - dX[1] / 2.0, nX[1] )
    X3_C = np.linspace( 0.0 + dX[2] / 2.0, 1.0 - dX[2] / 2.0, nX[2] )

    X1_F = np.linspace( 0.0, 1.2 * ( R + RT), nX[0] + 1 )
    X2_F = np.linspace( 0.0,     pi, nX[1] + 1 )
    X3_F = np.linspace( 0.0,    1.0, nX[2] + 1 )

    for j in range(nX[1]):
        for i in range(nX[0]):
            
            dis2 = R**2 + X1_C[i]**2 - 2.0 * R * X1_C[i] * m.sin(X2_C[j])
            #print(dis2)
            #print(RT**2)
            if (dis2 < RT**2 ):
                #print('Inside Ring')
                D[i,j] = D0
                V1_L[i,j] = -R0 * amp * W0 * m.sin( W0 * t )

    return D, V1_L, V2_L, V3_L, X1_C, X2_C, X3_C, \
           X1_F, X2_F, X3_F, dX


# Constants
cvel = 2.9979e10
Gc = 6.6726e-8

# --- Get user's HOME directory ---
HOME = subprocess.check_output( ["echo $HOME"], shell = True)
HOME = HOME[:-1].decode( "utf-8" ) + '/'

############# User input #############

#Root = '2D_3rd_896x128_COVn2_A2e-1_SW1e3_RIN1.16e3_Ent3.0'
Root = 'GR2D_M2.8_Mdot0.3_Rs120'
suffix = ''

DataDirectory = HOME + 'Research/' + '{:}/'.format( Root )
WritingDirectory = HOME + 'Research/Data/'
F = open( '{:}{:}{:}.txt'.format( WritingDirectory, Root, '_Data' ), 'w' )
PlotFileBaseName = '{:}{:}.plt'.format( Root, suffix )
print(PlotFileBaseName)
UsePhysicalUnits = True

CoordinateSystem = 'spherical'

#Analytical Test Inputs
AnalyticalStrain = False
if AnalyticalStrain:
    nX = [256, 256, 1]
    M0 = 1.989e34 #central mass
    ms = 1.989e33 #ring mass

    R0 = 6.0 * 2.0 * Gc * M0 / cvel**2 #equilibrium distance
    RT0 = 0.1 * R0 #thickness
    W0 = m.sqrt(Gc * M0 / R0**3) #angular frequency
    D0 = ms / (2.0 * pi**2 * RT0**2 * R0)

    amp = 0.2 #amplitude

# Observer angles

thobs = pio2
phobs = 0.0

# --- Get Data analytically for multiple time steps ---
if AnalyticalStrain == True:

    
    ktime = 200
    dt = 2.0 * pi * 2.0 / W0 / ktime
    Time = 0.0

    n2m_r = np.zeros( shape = ( ktime, 5, nX[0] ), dtype = np.cdouble  )

    Times = np.zeros( shape = (ktime) )
    
    Acoef = 16.0 * Gc * pi**1.5 / (cvel**4 * m.sqrt(15.0))
    A2 = np.zeros( shape = (ktime) )
    
    rh_const = Gc * ms * (R0 * W0)**2 * amp / cvel**4 
    rhA = np.zeros( shape = (ktime) )

    for iTime in range(ktime):

        D, V1_L, V2_L, V3_L, X1_C, X2_C, X3_C, X1_F, X2_F, X3_F, dX \
            = GetData_Analytic_Strain( cvel, Gc, pi, M0, ms, R0, RT0, W0, D0, amp, nX, Time )
        
        A2[iTime] = Acoef * ms * (R0 * W0)**2 * (amp * m.cos(W0 * Time) + amp**2 * m.cos(2.0 * W0 * Time)) / (2.0 * pi)

        rhA[iTime] = rh_const * ( m.cos(W0 * Time) + amp * m.cos(2.0 * W0 * Time) ) * m.sin(thobs)**2

        Time = Time + dt

        nX1 = int(nX[0])
        nX2 = int(nX[1])
        nX3 = int(nX[2])

        dX1 = dX[0]
        dX2 = dX[1]
        dX3 = dX[2]

        Times[iTime] = Time

        # Make one velocity array for passing into integmerid.

        V_L = np.zeros( shape = ( nX1, nX2, 3 ) )
        V_L[:,:,0] = V1_L
        V_L[:,:,1] = V2_L
        V_L[:,:,2] = V3_L
        
        #print(D)
        #print(V1_L)
        # Perform the azimuthal integration (trapezoid rule).
        if nX3 == 1: 
            X3_C = X3_C * twopi
            X3_F = X3_F * twopi
            dX3 = dX3 * twopi
    
        for k in range( 0, nX3 ): 
          
            #print(k)
            #print(X3_C[k], X3_F[k], X3_F[k+1], dX3)
            left  = np.zeros( shape = ( ktime, 5, nX[0] ), dtype = np.cdouble )
            right = np.zeros( shape = ( ktime, 5, nX[0] ), dtype = np.cdouble )
            
            # Perform the meridional integration.

            left[iTime,:,:] = meridional_int_module.integmerid( nx = nX1, ny = nX2, \
                                                                         x_cf = X1_C, x_ef = X1_F, \
                                                                         dx = dX1, \
                                                                         y_cf = X2_C, y_ef = X2_F, \
                                                                         dy = dX2, \
                                                                         radialzonespershell = 1, \
                                                                         numberofradialshells = nX1, \
                                                                         phi = X3_F[k], \
                                                                         rho0 = D, vel0 = V_L )

            right[iTime,:,:] = meridional_int_module.integmerid( nx = nX1, ny = nX2, \
                                                                         x_cf = X1_C, x_ef = X1_F, \
                                                                         dx = dX1, \
                                                                         y_cf = X2_C, y_ef = X2_F, \
                                                                         dy = dX2, \
                                                                         radialzonespershell = 1, \
                                                                         numberofradialshells = nX1, \
                                                                         phi = X3_F[k+1], \
                                                                         rho0 = D, vel0 = V_L )

            n2m_r[iTime,:,:] = ( dX3 / 2.0 ) * ( left[iTime,:,:] + right[iTime,:,:] )
            


    # Evaluate the amplitudes using central differencing.

    A2m_r   = np.zeros( shape = ( ktime, 5, nX[0] ), dtype = np.cdouble )
    A2m_tot = np.zeros( shape = ( ktime, 5 ), dtype = np.cdouble )

    for iTime in range( 1, ktime-1 ):

        dt = half * ( Times[iTime+1] - Times[iTime-1] )

        A2m_r[iTime,:,:] = half * ( n2m_r[iTime+1,:,:] - n2m_r[iTime-1,:,:] ) / dt

        A2m_tot[iTime,:] = np.sum( A2m_r[iTime,:,:], axis = 1 )
    print(A2)
    #plt.plot( Times[1:ktime-1], A2m_tot[1:ktime-1,2], 'x', color = 'black', markersize = 3 )
    #plt.plot( Times[:], np.real( A2m_tot[:,1] ), color = 'green' )
    #plt.plot( Times[:], np.real( A2m_tot[:,3] ), color = 'blue' )
    #plt.plot( Times[:], np.real( A2m_tot[:,4] ), color = 'orange' )
    #plt.plot( Times[:], np.real( A2m_tot[:,0] ), color = 'red' ) 

    #plt.plot( Times[:], A2, '-', color = 'red' )
    #plt.xlabel("Time(sec)")
    #plt.ylabel('$A_{20}$')
    #plt.legend( ['Numeric', 'Analytic'], loc = "upper left" )

    
    # Get the tensor spherical harmonics.

    f2m = np.zeros( shape = ( 5, 2 ), dtype = np.cdouble )

    f2m = basisf( thobs, phobs )

    # Calculate strains.

    rh_r   = np.zeros( shape = ( ktime, 2, nX[0] ), dtype = np.double )
    rh_tot = np.zeros( shape = ( ktime, 2 ), dtype = np.double )
    if nX3 == 1:
        for iTime in range( 1, ktime-1 ):
            print(iTime)
            
            rh_r[iTime,0,:] = rh_r[iTime,0,:] + np.real( A2m_r[iTime,2,:] * f2m[2,0] )
            rh_r[iTime,1,:] = rh_r[iTime,1,:] + np.real( A2m_r[iTime,2,:] * f2m[2,1] )

            rh_tot[iTime,0] = np.sum( rh_r[iTime,0,:], axis = 0 )
            rh_tot[iTime,1] = np.sum( rh_r[iTime,1,:], axis = 0 )
    

    else:
        for iTime in range( 1, ktime-1 ):
            print(iTime)

            for m in range( 0, 5 ):

                rh_r[iTime,0,:] = rh_r[iTime,0,:] + np.real( A2m_r[iTime,m,:] * f2m[m,0] )
                rh_r[iTime,1,:] = rh_r[iTime,1,:] + np.real( A2m_r[iTime,m,:] * f2m[m,1] )

            rh_tot[iTime,0] = np.sum( rh_r[iTime,0,:], axis = 0 )
            rh_tot[iTime,1] = np.sum( rh_r[iTime,1,:], axis = 0 )
    
    plt.plot( Times[1:ktime-2], rh_tot[1:ktime-2,0], color = 'black' )
    plt.plot( Times[:], rhA, color = 'red' )
    plt.xlabel("Postbounce Time (ms)")
    plt.ylabel("D h (cm)")
    plt.legend( ['Numeric', 'Analytic'], loc = "upper left" )
    #plt.ylim(-250, 250)

    #plt.plot( Times[1:ktime-2], rh_tot[1:ktime-2,0], color = 'black' )
    #plt.plot( Times[1:ktime-2], rh_tot[1:ktime-2,1], color = 'red' )
    #plt.xlabel("Postbounce Time (ms)")
    #plt.ylabel("D h (cm)")
    #plt.legend( ['$D h_x$','$D h_+$'], loc = "upper left" )
    #plt.plot( X1_C[:], rh_r[nFiles-2,0,:] )
    #plt.plot( X1_C[:], rh_r[nFiles-2,1,:] )
    plt.show()
# --- Get Data for files in an array ---
elif AnalyticalStrain == False:
    File, FileArray = ChoosePlotFile( DataDirectory, PlotFileBaseName, \
                                  ReturnFileArray = True, Verbose = False )

    D, V1_L, V2_L, V3_L, X1_C, X2_C, X3_C, X1_F, X2_F, X3_F, xL, xU, nX, dX, Time \
        = GetDataR( DataDirectory, PlotFileBaseName, File, FileArray, \
                    CoordinateSystem, UsePhysicalUnits, \
                    Verbose = False )

    iFile = 0
    nFiles = len(FileArray)

    n2m_r = np.zeros( shape = ( nFiles, 5, nX[0] ), dtype = np.cdouble  )

    Times = np.zeros( shape = ( nFiles ) )

    for File in FileArray:

        print(File)

        D, V1_L, V2_L, V3_L, X1_C, X2_C, X3_C, X1_F, X2_F, X3_F, xL, xU, nX, dX, Time \
            = GetDataR( DataDirectory, PlotFileBaseName, File, FileArray, \
                        CoordinateSystem, UsePhysicalUnits, \
                        Verbose = False )

        #print(D)
        print(Time)
        nX1 = int(nX[0])
        nX2 = int(nX[1])
        nX3 = int(nX[2])

        dX1 = dX[0]
        dX2 = dX[1]
        dX3 = dX[2]

        Times[iFile] = Time / 1000

        # Make one velocity array for passing into integmerid.

        V_L = np.zeros( shape = ( nX1, nX2, 3 ) )
        V_L[:,:,0] = V1_L
        V_L[:,:,1] = V2_L
        V_L[:,:,2] = V3_L

        # Perform the meridional integration.

        n2m_r[iFile,:,:] = twopi * meridional_int_module.integmerid( nx = nX1, ny = nX2, \
                                                                     x_cf = X1_C, x_ef = X1_F, \
                                                                     dx = dX1, \
                                                                     y_cf = X2_C, y_ef = X2_F, \
                                                                     dy = dX2, \
                                                                     radialzonespershell = 1, \
                                                                     numberofradialshells = nX1, \
                                                                     phi = pi, \
                                                                     rho0 = D, vel0 = V_L )

        iFile += 1

    # Evaluate the amplitudes using central differencing.

    A2m_r   = np.zeros( shape = ( nFiles, 5, nX[0] ), dtype = np.cdouble )
    A2m_tot = np.zeros( shape = ( nFiles, 5 ), dtype = np.cdouble )

    for iFile in range( 1, nFiles-1 ):
        print(iFile)

        dt = half * ( Times[iFile+1] - Times[iFile-1] )

        A2m_r[iFile,2,:] = half * ( n2m_r[iFile+1,2,:] - n2m_r[iFile-1,2,:] ) / dt

        A2m_tot[iFile,2] = np.sum( A2m_r[iFile,2,:], axis = 0 )

    # Get the tensor spherical harmonics.

    f2m = np.zeros( shape = ( 5, 2 ), dtype = np.cdouble )

    f2m = basisf( thobs, phobs )

    # Calculate strains.

    rh_r   = np.zeros( shape = ( nFiles, 2, nX[0] ), dtype = np.double )
    rh_tot = np.zeros( shape = ( nFiles, 2 ), dtype = np.double )
    
    for iFile in range( 1, nFiles-1 ):
        print(iFile)

        for m in range( 0, 5 ):

            rh_r[iFile,0,:] = rh_r[iFile,0,:] + np.real( A2m_r[iFile,m,:] * f2m[m,0] )
            rh_r[iFile,1,:] = rh_r[iFile,1,:] + np.real( A2m_r[iFile,m,:] * f2m[m,1] )

        rh_tot[iFile,0] = np.sum( rh_r[iFile,0,:], axis = 0 )
        rh_tot[iFile,1] = np.sum( rh_r[iFile,1,:], axis = 0 )

        F.write( str(Times[iFile]) + ' ' )
        F.write( str(rh_tot[iFile,0]) + ' '  )
        F.write( str(rh_tot[iFile,1]) + '\n' )

    F.close()
       
    plt.plot( Times[1:nFiles-2] * 1000, rh_tot[1:nFiles-2,0], color = 'black' )
    plt.plot( Times[1:nFiles-2] * 1000, rh_tot[1:nFiles-2,1], color = 'red' )
    plt.xlabel("Postbounce Time (ms)")
    plt.ylabel("D h (cm)")
    plt.legend( ['$D h_x$','$D h_+$'], loc = "upper left" )
    #plt.plot( X1_C[:], rh_r[nFiles-2,0,:] )
    #plt.plot( X1_C[:], rh_r[nFiles-2,1,:] )
    plt.show()
