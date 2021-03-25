#!/usr/local/bin/python3

import numpy as np
import os
import subprocess
import math
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, SymLogNorm

from ReadFieldsHDF import ReadFields

print( '' )
print( 'Running PlotFieldsHDF.py...' )
print( '---------------------------' )


############################ User Input ############################

# --- Define root path for data ---


suffix = './Output100/'
if(suffix=='./Output100/'):
    nT=300
if(suffix=='./Output200/'):
    nT=600
if(suffix=='./Output50/'):
    nT=150

FigTitle = suffix[:-1]



Field             = 'GF_233'
UseSemiLogYScale  = True

TimeUnit          = ''
LengthUnit        = ''
Dimension         = ['X1','X2']#,'X3']
UseSemiLogXScale  = False
Snapshots         = [0]
cmap              = 'Purples'
PlotTranspose     = False
PlotContours      = False
yScale            = 1.0
PlotTwoFields     = False

SaveInitFinal = False
SaveInitFinalName = 'InitFinal_{:}_{:}.dat'.format( suffix[18:], Field )

############################ End of User Input #####################

Snapshots = np.array( Snapshots )
Dimension = np.array( Dimension )

Polar            = True
CoordinateSystem = 'SPHERICAL'


PathToData =suffix

nFiles = Snapshots.shape[0]

names = ReadFields( PathToData, Snapshots )


##### Plotting #####

fig = plt.figure( figsize = (6,4) )
ax  = fig.add_subplot( 111, polar = Polar )
fig.suptitle( FigTitle )

if( Dimension.shape[0] == 1 ):

    X = names[Dimension[0]]

    # Which slice to plot
    iSS = 0

    alpha = np.linspace( 0.2, 1.0, Snapshots.shape[0] )

    if  ( Dimension[0] == 'X1' ):

        Data = np.empty( (Snapshots.shape[0],X.shape[0]), float )
        SqrtGm = np.empty( (Snapshots.shape[0],X.shape[0]), float )

        if( PlotTwoFields ):

            Data2 = np.empty( (Snapshots.shape[0],X.shape[0]), float )
            SqrtGm2 = np.empty( (Snapshots.shape[0],X.shape[0]), float )

        for i in range( nFiles ):

            Data  [i] = names[Field][i][0,:] / yScale
            ax.plot(X,Data[i])
            #SqrtGm[i] = names['GF_Sg'][i][0,iSS,:]


    elif( Dimension[0] == 'X2' ):

        Data = np.empty( (Snapshots.shape[0],X.shape[0]), float )

        for i in range( nFiles ):

            Data[i] = names[Field][i][0,:,iSS]

            ax.plot( X, Data[i], label = 't = {:.2e} {:}'.format \
              ( T[i], TimeUnit ) )

    if( SaveInitFinal ):

#        np.savetxt( SaveInitFinalName, \
#                    np.vstack( (X,SqrtGm[0],Data[0],SqrtGm[1],Data[1]) ), \
#                    fmt = '%.16e' )

        np.savetxt( 'Final_{:}_{:}.dat'.format( suffix[7:-1], Field ), \
                    Data[-1] )

    xlim = [ X[0], X[-1] ]
    ylim = [ Data.min(), Data.max() ]

    if( UseSemiLogXScale ): ax.set_xscale( 'log' )
    if( UseSemiLogYScale ): ax.set_yscale( 'log' )

    ax.set_xlim( xlim )
    ax.set_ylim( ylim )

    if( PlotTwoFields ):

        from matplotlib.lines import Line2D

        leg1 = ax.legend()
        custom_lines = [ Line2D( [0], [0], color = 'black' ), \
                         Line2D( [0], [0], color = 'red' ) ]
        ax.legend( custom_lines, [ '{:}'.format( Field ), \
                                   '{:}'.format( Field2 ) ] )


        ax.add_artist(leg1)

    else:

        ax.legend()
        ax.set_ylabel( Field )

    ax.set_xlabel( '{:} {:}'.format( Dimension[0], LengthUnit ) )

#    ax.axvline( 180.0, color = 'k' )
#    X = np.linspace( 40, 540, 129 )
#    for i in range( X.shape[0] ):
#        ax.axvline(X[i])

#elif( Dimension.shape[0] == 2 ):
else:

    X1 = names[Dimension[0]]
    X2 = names[Dimension[1]]

    Data   = np.empty( (Snapshots.shape[0],X2.shape[0],X1.shape[0]), float )
    SqrtGm = np.empty( (Snapshots.shape[0],X2.shape[0],X1.shape[0]), float )

    for i in range( nFiles ):

        if( PlotTranspose ):

            Data[i] \
              = np.abs( names[Field][i][0,:,:] - names[Field][i][0,:,:].T )

        else:

            Data[i]   = names[Field][i][:,:]
            #SqrtGm[i] = names['GF_Sg'][i][0,:,:]

    if( SaveInitFinal ):
        np.savetxt( 'Init_{:}_{:}.dat'.format( suffix[7:10], Field ),Data[0] )
        np.savetxt( 'Final_{:}_{:}.dat'.format( suffix[7:10], Field ),Data[-1] )

    # Which snapshot to plot
    iSS = -1

    vmin = Data[iSS].min()
    vmax = Data[iSS].max()

    Norm = None

    extent = [ X1[0], X1[-1], X2[0], X2[-1] ]

    if( UseSemiLogYScale ):

        Norm = LogNorm()

        if( vmin <= 0.0 ):

            Norm = SymLogNorm( linthresh = 1.0e1 )

    if( CoordinateSystem == 'SPHERICAL' ):

        """
        Taken from:
        https://brushingupscience.com/2016/06/21/
        matplotlib-animations-the-easy-way/
        """
        theta, r = np.meshgrid( X2, X1 )

        im = ax.pcolormesh( theta, r, Data[iSS,:,:].T, \
                            cmap = cmap, \
                            vmin = vmin, vmax = vmax, \
                            norm = Norm )

        ax.set_thetamin( 180.0/np.pi * X2[0 ] )
        ax.set_thetamax( 180.0/np.pi * X2[-1] )

        ax.set_theta_zero_location( 'W' ) # z-axis horizontal
        ax.set_theta_direction( -1 )



    else:

        if( PlotContours ):

            levels = np.linspace( vmin, vmax, 30 )

            im = ax.contour( Data[iSS], \
                             levels = levels, \
                             extent = extent, \
                             cmap   = cmap )

        else:

            im = ax.imshow( Data[iSS], \
                            vmin   = vmin, \
                            vmax   = vmax, \
                            norm   = Norm, \
                            extent = extent, \
                            origin = 'lower', \
                            cmap   = cmap )



    #cbar = fig.colorbar( im )
    #cbar.set_label( Field )

#plt.savefig( HOME + 'Desktop/{:}_{:}.png'.format( suffix[:-1], Field ) )
#plt.show()
#plt.close()




Analytic   = np.zeros([nT,300],dtype=np.float64)
Error   = np.zeros([300,nT],dtype=np.float64)

for i in range (0,300):
    for j in range(0,nT):
        if (Field=='GF_212' or Field=='GF_313'):
            Analytic[j][i]=1.0/X1[i]
        if(Field=='GF_122'):
            Analytic[j][i]=-X1[i]
        if(Field=='GF_133'):
            Analytic[j][i]=-X1[i]*math.sin(X2[j])*math.sin(X2[j])
        if(Field=='GF_233'):
            Analytic[j][i]=-math.sin(X2[j])*math.cos(X2[j])
        if(Field=='GF_323'):
            Analytic[j][i]=1.0/math.tan(X2[j])
        Error[i][j]=100.0*abs(Data[0][j][i]-Analytic[j][i])/abs(Analytic[j][i])



fig2=plt.figure(2)
a=ax.pcolormesh(X2,X1,Error,norm=LogNorm())
cbar = fig.colorbar( a )
cbar.set_label( Field )




"""
Root = '/Users/dunhamsj/Research/Codes/poseidon/Poseidon_Output/'

Results = Root + 'Results/'
Sources = Root + 'Objects/Sources/'


prR     = np.loadtxt( Results + 'Results_Radial_Locs_00001.out' ) / 1.0e5
pBeta1 = np.loadtxt( Results + 'Results_Beta1_00001.out' ) / 1.0e5
pPsi   = np.loadtxt( Results + 'Results_ConFactor_00001.out' )
pAlpha = np.loadtxt( Results + 'Results_Lapse_00001.out' )

print( 'SHAPE( Results ) = ', prR.shape )

prS = np.loadtxt( Sources + 'Sources_Radial_Locs_00001.out' ) / 1.0e5
pE  = np.loadtxt( Sources + 'Sources_E_00001.out' )
pS  = np.loadtxt( Sources + 'Sources_S_00001.out' )
pS1 = np.loadtxt( Sources + 'Sources_S1_00001.out' )

print( 'SHAPE( Sources ) = ', prS.shape )

pls = '.-'
pms = 5.0

fig = plt.figure()

ax1 = fig.add_subplot( 321 )
ax1.semilogx( prR, pBeta1, pls, markersize = pms, label = 'Poseidon' )
ax1.set_ylabel( r'$\beta^{1}$ [km/s]' )
ax1.get_xaxis().set_visible( False )
ax1.legend()

ax2 = fig.add_subplot( 323 )
ax2.semilogx( prR, pAlpha, pls, markersize = pms )
ax2.set_ylabel( r'$\alpha$' )
ax2.get_xaxis().set_visible( False )

ax3 = fig.add_subplot( 325 )
ax3.semilogx( prR, pPsi, pls, markersize = pms )
ax3.set_ylabel( r'$\psi$' )
ax3.set_xlabel( 'Radial Coordinate [km]' )

ax4 = fig.add_subplot( 322 )
ax4.semilogx( prS, pE, pls, markersize = pms )
ax4.set_ylabel( r'$E$' )
ax4.get_xaxis().set_visible( False )

ax5 = fig.add_subplot( 324 )
ax5.semilogx( prS, pS, pls, markersize = pms )
ax5.set_ylabel( r'$S$' )
ax5.get_xaxis().set_visible( False )

ax6 = fig.add_subplot( 326 )
ax6.semilogx( prS, pS1, pls, markersize = pms )
ax6.set_ylabel( r'$S^{1}$' )
ax6.set_xlabel( 'Radial Coordinate [km]' )

plt.subplots_adjust( hspace = 0.0 )
plt.show()

"""

os.system( 'rm -rf __pycache__' )
