from CellAverages import CellAverage
import rcparams
import h5py as h5
import matplotlib.pyplot as plt
import os
import glob
import sys
import matplotlib.colors as colors
import numpy as np
from matplotlib.ticker import AutoMinorLocator, LogLocator, LogFormatterMathtext

DirName = sys.argv[1]
DirName = DirName.removeprefix("../../../SandBox/dgExperiments_MHD_Relativistic_IDEAL/Output/")

DataFiles = sorted(glob.glob(sys.argv[1] + '/ShearingDisk_MagnetofluidFields*.h5'))
GeometryFile = sys.argv[1] + '/ShearingDisk_GeometryFields_000000.h5'

gf = h5.File(GeometryFile, 'r')

X1 = gf['Spatial Grid']['X1'][:] * 1.0e5
X2 = gf['Spatial Grid']['X2'][:] * 1.0e5
X3 = gf['Spatial Grid']['X3'][:]

X1_C = gf['Spatial Grid']['X1_C'][:] * 1.0e5
X2_C = gf['Spatial Grid']['X2_C'][:] * 1.0e5
X3_C = gf['Spatial Grid']['X3_C'][:]

dX1  = gf['Spatial Grid']['dX1'][:] * 1.0e5
dX2  = gf['Spatial Grid']['dX2'][:] * 1.0e5
dX3  = gf['Spatial Grid']['dX3'][:]

if( np.size(dX1) > 1):
	nDim = 1
	if( np.size(dX2) > 1):
		nDim = 2
		if( np.size(dX3) > 1):
			nDim = 3

nNodes = int(np.size(X1)) // int(np.size(X1_C))

Psi    = gf['Geometry Fields']['Conformal Factor'][:,:,:]
Psi4   = Psi**4
Psi6   = Psi**6

Gm11   = CellAverage( nDim, nNodes, gf['Geometry Fields']['Spatial Metric Component (11)'][:,:,:] )
Gm22   = CellAverage( nDim, nNodes, gf['Geometry Fields']['Spatial Metric Component (22)'][:,:,:] )
Gm33   = CellAverage( nDim, nNodes, gf['Geometry Fields']['Spatial Metric Component (33)'][:,:,:] * 1.0e10 )
SqrtGm = CellAverage( nDim, nNodes, gf['Geometry Fields']['Sqrt Spatial Metric Determinant'][:,:,:] * 1.0e5 )

nFiles = len(DataFiles)

Means = [np.zeros(nFiles), np.zeros(nFiles), np.zeros(nFiles), np.zeros(nFiles), np.zeros(nFiles), np.zeros(nFiles), np.zeros(nFiles)]

def int_total_1d(quant):
	int_total = 2.0 * np.pi * 1.0e5 * np.sum( dX1[:] * quant[0,0,:] * SqrtGm[0,0,:] )
	return int_total

Vol = int_total_1d( np.ones_like( SqrtGm[:,:,:] ) )
VolChk = np.pi * 1.0e5 * (X1_C[-1]**2 - X1_C[0]**2)

Times = np.zeros(nFiles)

iFile = 0
for file in DataFiles:

	f = h5.File(file, 'r')

	B1 = CellAverage( nDim, nNodes, f['Magnetofluid Fields']['Conserved']['Conserved Magnetic Field (1)'][:,:,:])
	B2 = CellAverage( nDim, nNodes, f['Magnetofluid Fields']['Conserved']['Conserved Magnetic Field (2)'][:,:,:])
	B3 = CellAverage( nDim, nNodes, f['Magnetofluid Fields']['Conserved']['Conserved Magnetic Field (3)'][:,:,:])
	IntE = CellAverage( nDim, nNodes, f['Magnetofluid Fields']['Primitive']['Internal Energy Density'][:,:,:])

	D  = CellAverage( nDim, nNodes, f['Magnetofluid Fields']['Primitive']['Comoving Baryon Density'][:,:,:])
	V1 = CellAverage( nDim, nNodes, f['Magnetofluid Fields']['Primitive']['Three-Velocity (1)'][:,:,:] * 1.0e5)
	V2 = CellAverage( nDim, nNodes, f['Magnetofluid Fields']['Primitive']['Three-Velocity (2)'][:,:,:] * 1.0e5)
	V3 = CellAverage( nDim, nNodes, f['Magnetofluid Fields']['Primitive']['Three-Velocity (3)'][:,:,:])

	B1Sq = Gm11[:,:,:] * B1[:,:,:]**2
	B2Sq = Gm22[:,:,:] * B2[:,:,:]**2
	B3Sq = Gm33[:,:,:] * B3[:,:,:]**2
	KE1  = D[:,:,:] * Gm11[:,:,:] * V1[:,:,:]**2 / 2.0
	KE2  = D[:,:,:] * Gm22[:,:,:] * V2[:,:,:]**2 / 2.0
	KE3  = D[:,:,:] * Gm33[:,:,:] * V3[:,:,:]**2 / 2.0

	Means[0][iFile] = int_total_1d( B1Sq ) / ( 2.0 * Vol )
	Means[1][iFile] = int_total_1d( B2Sq ) / ( 2.0 * Vol )
	Means[2][iFile] = int_total_1d( B3Sq ) / ( 2.0 * Vol )
	Means[3][iFile] = int_total_1d( IntE ) / Vol
	Means[4][iFile] = int_total_1d( KE1 ) / Vol
	Means[5][iFile] = int_total_1d( KE2 ) / Vol
	Means[6][iFile] = int_total_1d( KE3 ) / Vol

	Times[iFile] = f['Time'][0]

	print( 'File ' + str(iFile) + ' of ' + str(nFiles) )

	iFile = iFile + 1

fig, ax = plt.subplots(1, 1, figsize=(10,7))

ax.semilogy( Times[:], Means[0][:], 'r--', linewidth = 3.25, label = r'$\mathrm{B_r}$' )
ax.semilogy( Times[:], Means[1][:], 'b--', linewidth = 3.25, label = r'$\mathrm{B_z}$')
ax.semilogy( Times[:], Means[2][:], 'g--', linewidth = 3.25, label = r'$\mathrm{B_\theta}$' )
ax.semilogy( Times[:], Means[3][:], 'm-', linewidth = 3.25, label = 'Thermal')
ax.semilogy( Times[:], Means[4][:], 'r:', linewidth = 3.25, label = r'$\mathrm{KE_r}$' )
ax.semilogy( Times[:], Means[5][:], 'b:', linewidth = 3.25, label = r'$\mathrm{KE_z}$')
ax.semilogy( Times[:], Means[6][:], 'g:', linewidth = 3.25, label = r'$\mathrm{KE_\theta}$' )

ax.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,), numticks = 100))
ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2,10), numticks = 100))
ax.yaxis.set_major_formatter(LogFormatterMathtext())
ax.xaxis.set_minor_locator(AutoMinorLocator(5))

ax.set_title( DirName, fontsize = 12 )
ax.set_ylabel( r'Mean Energy Density [$\mathrm{erg ~ cm^{-3}}$]' )
ax.set_xlabel( r'Time [$\mathrm{ms}$]' )
ax.set_ylim( 1.0e25, 1.0e35 )
ax.set_xlim( 0.0, Times[-1] )
ax.legend(bbox_to_anchor=(1.05,1), loc='upper left')
fig.tight_layout()
fig.savefig( DirName + '_Energetics.png', dpi = 250 )
