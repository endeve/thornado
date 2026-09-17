from mpi4py import MPI
import gc
import numpy as np
import yt
from os.path import isfile
import matplotlib.pyplot as plt
from sys import argv
import warnings
import scipy.integrate as integrate
import sys
sys.path.append('../../AMReX')
import piecewise_regression
from Utilities.Files import GetFileNumberArray

yt.set_log_level(50)
yt.enable_parallelism()

PlotDirectory = argv[1] + '/'
PlotBaseName = 'ShearingDisk2D.plt.'

FileNumberArray \
= GetFileNumberArray \
  ( PlotDirectory, \
    PlotBaseName, \
    -1, -1, 1 )

#size = 12
size = np.size(FileNumberArray)

FileNames = []

for index in range(0, size, 1):

    FileNames.append( PlotDirectory + '{:}{:}'.format( PlotBaseName, \
                      str( FileNumberArray[index] ).zfill( 8 ) ) )

# Let yt build timeseries for parallel processing.

ts = yt.DatasetSeries(FileNames)

data = {}

for sto, ds in ts.piter(storage=data):

    ad = ds.all_data()

    t = ds.current_time.to_ndarray()

    print('Time: ' + str(t))

    nX1 = ds.domain_dimensions[0]
    nX2 = ds.domain_dimensions[1]
    nX3 = ds.domain_dimensions[2]

    dX1   = ad[('boxlib', 'dX1')].to_ndarray().reshape( nX1, nX2, nX3 )[:,0,0] * 1.0e5
    dX2   = ad[('boxlib', 'dX2')].to_ndarray().reshape( nX1, nX2, nX3 )[0,:,0] * 1.0e5
    dX3   = ad[('boxlib', 'dX3')].to_ndarray().reshape( nX1, nX2, nX3 )[0,0,:]

    Gmdd11 = ad[('boxlib', 'GF_Gm_11' )].to_ndarray().reshape( nX1, nX2, nX3 )
    Gmdd22 = ad[('boxlib', 'GF_Gm_22' )].to_ndarray().reshape( nX1, nX2, nX3 )
    Gmdd33 = ad[('boxlib', 'GF_Gm_33' )].to_ndarray().reshape( nX1, nX2, nX3 ) * 1.0e10
    SqrtGm = ad[('boxlib', 'GF_SqrtGm')].to_ndarray().reshape( nX1, nX2, nX3 ) * 1.0e5

    Vol = np.sum( np.sum( np.sum( SqrtGm[:,:,:] * dX3[:], axis = 2 ) * dX2[:], axis = 1 ) * dX1[:], axis = 0 )

    PD  = ad[('boxlib', 'PM_D' )].to_ndarray().reshape( nX1, nX2, nX3 )
    PB1 = ad[('boxlib', 'PM_B1')].to_ndarray().reshape( nX1, nX2, nX3 )
    PB2 = ad[('boxlib', 'PM_B2')].to_ndarray().reshape( nX1, nX2, nX3 )
    PB3 = ad[('boxlib', 'PM_B3')].to_ndarray().reshape( nX1, nX2, nX3 )
    CB1 = ad[('boxlib', 'CM_B1')].to_ndarray().reshape( nX1, nX2, nX3 )
    CB2 = ad[('boxlib', 'CM_B2')].to_ndarray().reshape( nX1, nX2, nX3 )
    CB3 = ad[('boxlib', 'CM_B3')].to_ndarray().reshape( nX1, nX2, nX3 )
    PV1 = ad[('boxlib', 'PM_V1')].to_ndarray().reshape( nX1, nX2, nX3 ) * 1.0e5
    PV2 = ad[('boxlib', 'PM_V2')].to_ndarray().reshape( nX1, nX2, nX3 ) * 1.0e5
    PV3 = ad[('boxlib', 'PM_V3')].to_ndarray().reshape( nX1, nX2, nX3 )
    PE  = ad[('boxlib', 'PM_E') ].to_ndarray().reshape( nX1, nX2, nX3 )

    PB1E = Gmdd11 * PB1**2 / 2.0
    PB2E = Gmdd22 * PB2**2 / 2.0
    PB3E = Gmdd33 * PB3**2 / 2.0
    CB1E = Gmdd11 * CB1**2 / 2.0
    CB2E = Gmdd22 * CB2**2 / 2.0
    CB3E = Gmdd33 * CB3**2 / 2.0
    PV1E = ( 1.0 / 2.0 ) * PD * Gmdd11 * PV1**2
    PV2E = ( 1.0 / 2.0 ) * PD * Gmdd22 * PV2**2
    PV3E = ( 1.0 / 2.0 ) * PD * Gmdd33 * PV3**2

    AvgE = np.zeros( 10 )

    AvgE[0] = ( 1.0 / Vol ) * np.sum( np.sum( np.sum( SqrtGm[:,:,:] * PB1E[:,:,:] * dX3[:], axis = 2 ) * dX2[:], axis = 1 ) * dX1[:], axis = 0 )
    AvgE[1] = ( 1.0 / Vol ) * np.sum( np.sum( np.sum( SqrtGm[:,:,:] * PB2E[:,:,:] * dX3[:], axis = 2 ) * dX2[:], axis = 1 ) * dX1[:], axis = 0 )
    AvgE[2] = ( 1.0 / Vol ) * np.sum( np.sum( np.sum( SqrtGm[:,:,:] * PB3E[:,:,:] * dX3[:], axis = 2 ) * dX2[:], axis = 1 ) * dX1[:], axis = 0 )
    AvgE[3] = ( 1.0 / Vol ) * np.sum( np.sum( np.sum( SqrtGm[:,:,:] * CB1E[:,:,:] * dX3[:], axis = 2 ) * dX2[:], axis = 1 ) * dX1[:], axis = 0 )
    AvgE[4] = ( 1.0 / Vol ) * np.sum( np.sum( np.sum( SqrtGm[:,:,:] * CB2E[:,:,:] * dX3[:], axis = 2 ) * dX2[:], axis = 1 ) * dX1[:], axis = 0 )
    AvgE[5] = ( 1.0 / Vol ) * np.sum( np.sum( np.sum( SqrtGm[:,:,:] * CB3E[:,:,:] * dX3[:], axis = 2 ) * dX2[:], axis = 1 ) * dX1[:], axis = 0 )
    AvgE[6] = ( 1.0 / Vol ) * np.sum( np.sum( np.sum( SqrtGm[:,:,:] * PV1E[:,:,:] * dX3[:], axis = 2 ) * dX2[:], axis = 1 ) * dX1[:], axis = 0 )
    AvgE[7] = ( 1.0 / Vol ) * np.sum( np.sum( np.sum( SqrtGm[:,:,:] * PV2E[:,:,:] * dX3[:], axis = 2 ) * dX2[:], axis = 1 ) * dX1[:], axis = 0 )
    AvgE[8] = ( 1.0 / Vol ) * np.sum( np.sum( np.sum( SqrtGm[:,:,:] * PV3E[:,:,:] * dX3[:], axis = 2 ) * dX2[:], axis = 1 ) * dX1[:], axis = 0 )
    AvgE[9] = ( 1.0 / Vol ) * np.sum( np.sum( np.sum( SqrtGm[:,:,:] * PE[:,:,:]   * dX3[:], axis = 2 ) * dX2[:], axis = 1 ) * dX1[:], axis = 0 )

    sto.result    = (t, AvgE)

    del ds

    gc.collect()

if yt.is_root():

    sl = np.empty(6)
    ic = np.empty(6)

    dlist = list(data.values())

    Times = np.array([tup[0] for tup in dlist])
    Avgs  = np.array([tup[1] for tup in dlist])

    # Fitting machinery (find way to do lims automatically?).

    FitLimsRad = [ 1.0e+22, 1.0e26 ]
    FitLimsRot = [ 1.0e+24, 5.0e27 ]

    RadFit = False
    RotFit = False

    IndexRad = np.nonzero( ( Avgs[:,0] > FitLimsRad[0] ) & ( Avgs[:,0] < FitLimsRad[1] ) )
    IndexRot = np.nonzero( ( Avgs[:,2] > FitLimsRot[0] ) & ( Avgs[:,2] < FitLimsRot[1] ) )

    pwRadFit = True
    pwRotFit = True

    if( np.size(IndexRad) > 0 ):
        RadFit = True
        FitTimes = Times[IndexRad[:]]
        FitAvgs  = Avgs[IndexRad[:], 0]
        sl[0], ic[0] = np.polyfit( FitTimes[:], np.log(FitAvgs[0]), 1)
    else:
        RadFit = False

    if( np.size(IndexRot) > 0 ):
        RotFit = True
        FitTimes = Times[IndexRot[:]]
        FitAvgs  = Avgs[IndexRot[:], 2]
        sl[2], ic[2] = np.polyfit( FitTimes[:], np.log(FitAvgs[0]), 1)
    else:
        RotFit = False

    # Sorting to avoid x being non-monotonic when running in parallel.

    sort_ind       = np.argsort(Times)
    Times          = Times[sort_ind]

    # Seperate attempt with piecewise regression library (https://joss.theoj.org/papers/10.21105/joss.03859)

    if( pwRadFit ):

        RadFit_pw     = piecewise_regression.Fit( Times[:], np.log(Avgs[sort_ind,0]), n_breakpoints = 2 )
        RadFit_pw_out = RadFit_pw.get_results()
        RadFit_pw.summary()
        RadFit_pw_const = RadFit_pw_out["estimates"]["const"]["estimate"] + RadFit_pw_out["estimates"]["alpha1"]["estimate"] * RadFit_pw_out["estimates"]["breakpoint1"]["estimate"]
        RadFit_pw_slope = RadFit_pw_out["estimates"]["alpha2"]["estimate"]

    if( pwRotFit ):

        RotFit_pw     = piecewise_regression.Fit( Times[:], np.log(Avgs[sort_ind,2]), n_breakpoints = 2 )
        RotFit_pw_out = RotFit_pw.get_results()
        RotFit_pw.summary()
        RotFit_pw_const = RotFit_pw_out["estimates"]["const"]["estimate"] + RotFit_pw_out["estimates"]["alpha1"]["estimate"] * RotFit_pw_out["estimates"]["breakpoint1"]["estimate"]
        RotFit_pw_slope = RotFit_pw_out["estimates"]["alpha2"]["estimate"]

    for i in range( 0, 10, 1 ):
        Avgs[:,i] = Avgs[sort_ind,i]

    E = np.hstack( (Times[:, np.newaxis], Avgs) )

    np.savetxt  ( PlotDirectory + 'Energetics.txt', E )

    plt.xlim( 0.0, 10.0 )
    plt.ylim( 1.0e+16, 10.0 * np.max(Avgs) )
    plt.xlabel('Time')
    plt.ylabel('Mean Energy Density ' + r'[$\mathrm{erg} ~ \mathrm{cm}^{-3}$]')

    plt.semilogy( Times[:], Avgs[:,0], color = 'red',    linestyle = '--', label = 'P. Rad.  Mag.')
    plt.semilogy( Times[:], Avgs[:,1], color = 'green',  linestyle = '--', label = 'P. Vert. Mag.')
    plt.semilogy( Times[:], Avgs[:,2], color = 'blue',   linestyle = '--', label = 'P. Rot.  Mag.')
    plt.semilogy( Times[:], Avgs[:,6], color = 'red',    linestyle = ':' , label = 'Rad.  KE'     )
    plt.semilogy( Times[:], Avgs[:,7], color = 'green',  linestyle = ':' , label = 'Vert. KE'     )
    plt.semilogy( Times[:], Avgs[:,8], color = 'blue',   linestyle = ':' , label = 'Rot.  KE'     )
    plt.semilogy( Times[:], Avgs[:,9], color = 'orange', linestyle = '-.', label = 'Therm.'       )

    if( False ):
         plt.plot( Times[:], np.exp( ic[0] ) * np.exp( sl[0] * Times[:] ), color = 'red',  linestyle = '-', label = 'Rad. Mag. Fit - sl: ' + str(sl[0]) )
    if( False ):
         plt.plot( Times[:], np.exp( ic[2] ) * np.exp( sl[2] * Times[:] ), color = 'blue', linestyle = '-', label = 'Rot. Mag. Fit - sl: ' + str(sl[2]) )

    if( pwRadFit ):
        plt.plot( Times[:], np.exp( RadFit_pw_const ) * np.exp( RadFit_pw_slope * Times[:] ), color = 'red',  linestyle = '-', label = 'Rad. Mag. Fit - sl: ' + str(RadFit_pw_slope) )
        RadFit_pw.plot_breakpoint_confidence_intervals(color = 'red')
    if( pwRotFit ):
        plt.plot( Times[:], np.exp( RotFit_pw_const ) * np.exp( RotFit_pw_slope * Times[:] ), color = 'blue', linestyle = '-', label = 'Rot. Mag. Fit - sl: ' + str(RotFit_pw_slope) )
        RotFit_pw.plot_breakpoint_confidence_intervals(color = 'blue')

    plt.legend(bbox_to_anchor=(1.05, 0.5), loc='center left')
    plt.savefig(PlotDirectory + 'Energetics.png', bbox_inches='tight')
