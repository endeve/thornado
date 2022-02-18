#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

def GetData( nN, nX, Grid, Reflux = '' ):

    i = np.loadtxt( 'nN{:}_nX{:}_{:}_I{:}.dat'.format( nN, nX, Grid, Reflux ) )
    f = np.loadtxt( 'nN{:}_nX{:}_{:}_F{:}.dat'.format( nN, nX, Grid, Reflux ) )

    N = np.int64( nN )

    iLevel = i[:,0]
    iX1    = i[:,1]
    iX2    = i[:,2]
    iX3    = i[:,3]
    dx1i   = i[:,4]
    dx2    = i[:,5]
    dx3    = i[:,6]
    x1i    = i[:,7:7+N].flatten()
    x2i    = i[:,7+N].flatten()
    x3i    = i[:,7+N+1].flatten()
    ui     = i[:,7+N+2:].flatten()
    di     = np.vstack( (x1i,ui) )
    indi   = di[0].argsort() # sort by x-value

    iLevel = f[:,0]
    iX1    = f[:,1]
    iX2    = f[:,2]
    iX3    = f[:,3]
    dx1f   = f[:,4]
    dx2    = f[:,5]
    dx3    = f[:,6]
    x1f    = f[:,7:7+N].flatten()
    x2f    = f[:,7+N].flatten()
    x3f    = f[:,7+N+1].flatten()
    uf     = f[:,7+N+2:].flatten()
    df     = np.vstack( (x1f,uf) )
    indf   = df[0].argsort() # sort by x-value

    return dx1i, dx1f, di, df, indi, indf

def PlotData \
    ( axs, di, df, indi, indf, c = 'k', label = '' ):

    axs[0].plot( di[0,indi], di[1,indi], c+'o', label = label + '_i' )
    axs[0].plot( df[0,indf], df[1,indf], c+'s', label = label + '_f' )
    axs[1].plot( di[0,indi], \
                  np.abs( ( df[1,indf] - di[1,indi] ) / di[1,indi] ), \
                  c, label = label )

    return

fig, axs = plt.subplots(2,1)

#dxiU, dxfU, diU, dfU, indiU, indfU \
#  = GetData( '02', '032', 'UniGrid' )
#PlotData( axs, diU, dfU, indiU, indfU, 'r', 'UniGrid' )

dxiU, dxfU, diU, dfU, indiU, indfU \
  = GetData( '02', '032', 'MultiGrid', Reflux = '_RefluxOn' )
PlotData( axs, diU, dfU, indiU, indfU, 'm', 'MultiGrid_RefluxOn' )

#dxiU, dxfU, diU, dfU, indiU, indfU \
#  = GetData( '02', '032', 'MultiGrid', Reflux = '_RefluxOff' )
#PlotData( axs, diU, dfU, indiU, indfU, 'b', 'MultiGrid_RefluxOff' )

axs[0].legend()
axs[1].legend()
plt.show()
plt.close()
