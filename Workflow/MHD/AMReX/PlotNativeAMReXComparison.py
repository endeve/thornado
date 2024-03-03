#!/usr/bin/env python3

import numpy as np
import yt
import h5py
import matplotlib.pyplot as plt

from Utilities.Files import GetFileNumberArray

AMReXPlotDirectory = '/home/jbuffal/thornado_MHD_3D/SandBox/AMReX/MHD_Relativistic_IDEAL/'

NativePlotDirectory = '/home/jbuffal/thornado_MHD_3D/SandBox/dgExperiments_MHD_Relativistic_IDEAL/Output/'

AMReXPlotFile  = 'plt_00006401'
NativePlotFile = 'OrszagTang2D_MagnetofluidFields_000101.h5'

NativeFieldType = 'Conserved'
NativeField     = 'Conserved Magnetic Field (1)'

AMReXField      = 'CM_B1'

ds = yt.load( AMReXPlotDirectory + AMReXPlotFile )
ds.force_periodicity()

NativeFile = h5py.File( NativePlotDirectory + NativePlotFile, 'r')

print(NativeFile['Magnetofluid Fields']['Primitive'].keys())
print(NativeFile['Magnetofluid Fields']['Conserved'].keys())
print(NativeFile['Magnetofluid Fields']['Diagnostic'].keys())

NativeData = np.transpose(NativeFile['Magnetofluid Fields'][NativeFieldType][NativeField][0,:,:])

MaxLevel = ds.index.max_level
nX       = ds.domain_dimensions
xL       = ds.domain_left_edge
xH       = ds.domain_right_edge

if ( nX[1] == 1 and nX[2] == 1 ):
    nDims = 1
elif( nX[2] == 1 ):
    nDims = 2
else:
    nDims = 3

CoveringGrid = ds.covering_grid( level = MaxLevel, left_edge = xL, dims = nX * 2**MaxLevel, num_ghost_zones = nX[0] )

AMReXData = CoveringGrid[AMReXField].to_ndarray()[:,:,0]

print('Comparing ' + AMReXField)
print( np.min( np.abs( NativeData - AMReXData ) ) )
print( np.max( np.abs( NativeData - AMReXData ) ) )

plt.figure(1)
plt.imshow(NativeData)
plt.figure(2)
plt.imshow(AMReXData)
plt.figure(3)
plt.imshow(NativeData - AMReXData)
plt.colorbar()
plt.clim(-1.0e-15, 1.0e-15)
plt.show()
