#!/usr/bin/env python3

import sys
import numpy as np
import yt
import h5py
import matplotlib.pyplot as plt

from matplotlib import animation
from Utilities.Files import GetFileNumberArray

sys.path.append('../')

PlotDirectory = '/home/jbuffal/thornado_MHD_3D/SandBox/AMReX/MHD_Relativistic_IDEAL/OT2D_2ndOrder_n128x128x1_SLimiting_BTVD1.75e00_CleaningOff_Damp0.00e00_PowellSource_Optimize/'

PlotBaseName      = 'plt_'

FileArray = GetFileNumberArray( PlotDirectory, PlotBaseName, -1, -1, 1 )

NumFiles = len(FileArray)

TimeVec = np.zeros(NumFiles)
TotalDivVec = np.zeros(NumFiles)

for i, FileNumber in enumerate(FileArray):

	PlotFile = PlotBaseName + '{:}'.format( str(FileNumber).zfill(8) )

	ds = yt.load( PlotDirectory + PlotFile )
	ds.force_periodicity()

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
	Time = ds.current_time.to_ndarray()

	PD   = CoveringGrid['PM_D'  ].to_ndarray()[:,:,:]
	PV1  = CoveringGrid['PM_V1' ].to_ndarray()[:,:,:]
	PV2  = CoveringGrid['PM_V2' ].to_ndarray()[:,:,:]
	PV3  = CoveringGrid['PM_V3' ].to_ndarray()[:,:,:]
	PE   = CoveringGrid['PM_E'  ].to_ndarray()[:,:,:]
	PB1  = CoveringGrid['PM_B1' ].to_ndarray()[:,:,:]
	PB2  = CoveringGrid['PM_B2' ].to_ndarray()[:,:,:]
	PB3  = CoveringGrid['PM_B3' ].to_ndarray()[:,:,:]
	CB1  = CoveringGrid['CM_B1' ].to_ndarray()[:,:,:]
	CB2  = CoveringGrid['CM_B2' ].to_ndarray()[:,:,:]
	CB3  = CoveringGrid['CM_B3' ].to_ndarray()[:,:,:]
	DDiv = CoveringGrid['DM_Div'].to_ndarray()[:,:,:]
	AGm  = CoveringGrid['AM_Gm' ].to_ndarray()[:,:,:]

	VSq = PV1 * PV1 + PV2 * PV2 + PV3 * PV3
	W = 1.0 / np.sqrt( 1.0 - VSq )
	bu0 = W * ( CB1 * PV1 + CB2 * PV2 + CB3 * PV3 )
	bd0 = -bu0

	PBSqr = bu0 * bd0 + PB1 * PB1 + PB2 * PB2 + PB3 * PB3

	MagEnergy   = (W * W - 0.5) * PBSqr - bu0 * bu0
	Div         = DDiv
	ThermEnergy = PE
	KinEnergy   = PD * W * W + AGm * PE * W * W
	AllEnergy   = MagEnergy + ThermEnergy + KinEnergy

	TotalMagEnergy   = np.sum(MagEnergy)
	TotalDiv         = np.sum(abs(Div))
	TotalThermEnergy = np.sum(ThermEnergy)
	TotalKinEnergy   = np.sum(KinEnergy)
	TotalAllEnergy   = np.sum(AllEnergy)

	TotalDivVec[i]      = TotalDiv
	TimeVec[i]          = Time

plt.plot(TimeVec, TotalDivVec)
plt.show()
