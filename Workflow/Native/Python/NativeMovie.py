from CellAverages import CellAverage
import rcparams
import h5py as h5
import matplotlib.pyplot as plt
import os
import glob
import sys
import matplotlib.colors as colors
import numpy as np
from matplotlib.ticker import SymmetricalLogLocator, AutoMinorLocator, LogFormatterMathtext

UseSymLog = False
Units = ''
lims = []
linthresh = 1.0e0

if sys.argv[3] == 'Conserved Magnetic Field (1)':
	FieldName = 'CM_B1'
	lims = [-1.0e15, 1.0e15]
	UseSymLog = True
	linthresh = 1.0e+10
	Units = r'CM_B1 $[\mathrm{G}]$'
elif sys.argv[3] == 'Conserved Magnetic Field (2)':
	FieldName = 'CM_B2'
	lims = [-1.0e15, 1.0e15]
	UseSymLog = True
	linthresh = 1.0e+10
	Units = r'CM_B2 $[\mathrm{G}]$'
elif sys.argv[3] == 'Conserved Magnetic Field (3)':
	FieldName = 'CM_B3'
	lims = [-1.0e9, 1.0e9]
	UseSymLog = True
	linthresh = 1.0e+04
	Units = r'CM_B3 $[\mathrm{G} ~ \mathrm{cm}^{-1}]$'
elif sys.argv[3] == 'gammaMRI':
	FieldName = 'gammaMRI'
elif sys.argv[3] == 'Three-Velocity (1)':
	FieldName = 'PM_V1'
	lims = [ -1.0e5, 1.0e5 ]
	UseSymLog = True
	Units = r'PM_V1 $[\mathrm{km} ~ \mathrm{s}^{-1}]$'
elif sys.argv[3] == 'Three-Velocity (3)':
	FieldName = 'PM_V3'
	lims = [ 1000, 3000 ]
	Units = r'PM_V3 $[\mathrm{s}^{-1}]$'
elif sys.argv[3] == 'Comoving Baryon Density':
	FieldName = 'PM_D'
	lims = [ 1.0e+13, 5.0e+13 ]
	Units = r'PM_D $[\mathrm{g} ~ \mathrm{cm}^{-3}]$'
elif sys.argv[3] == 'Conserved Momentum Density (1)':
	FieldName = 'CM_S1'
elif sys.argv[3] == 'Conserved Momentum Density (3)':
	FieldName = 'CM_S3'
	lims = [ 1.0e+30, 1.0e+32 ]
DirName = sys.argv[1]

DirName = DirName.removeprefix("../../../SandBox/dgExperiments_MHD_Relativistic_IDEAL/Output/")

DirPath = os.getcwd() + '/Images_' + DirName + '_' + FieldName

if not os.path.exists( DirPath ):
	os.mkdir( DirPath )
else:
	os.system( "rm -r " + DirPath )
	os.mkdir( DirPath )

DataFiles = sorted(glob.glob(sys.argv[1] + '/*MagnetofluidFields*.h5'))

iFile = 0
for file in DataFiles:

	f = h5.File(file, 'r')
	if sys.argv[2] == 'Conserved' or sys.argv[2] == 'Primitive' or sys.argv[2] == 'Auxiliary' or sys.argv[2] == 'Diagnostic':
		dset = f['Magnetofluid Fields'][sys.argv[2]][sys.argv[3]]
	elif sys.argv[2] == 'Post-Processed' and sys.argv[3] == 'gammaMRI':
		dset = f['Magnetofluid Fields']['Primitive']['Three-Velocity (3)']

	X1 = f['Spatial Grid']['X1']
	X2 = f['Spatial Grid']['X2']
	X3 = f['Spatial Grid']['X3']

	time = f['Time']

	if dset.shape[2] > 1:
		dim = 1
		if dset.shape[1] > 1:
			dim = 2
			if dset.shape[0] > 1:
				dim = 3

	if dim == 1:

		fig, ax = plt.subplots(1, 1, figsize=(10,7))
		if sys.argv[2] == 'Post-Processed':
			if sys.argv[3] == 'gammaMRI':
				q = - ( X1[:] / dset[0,0,:] ) * np.gradient( dset[0,0,:], X1[:], axis = 0, edge_order = 2 )
				gammaMRI = q[:] * dset[0,0,:]
				ax.plot(X1[:], gammaMRI[:])
				ax.set_yscale('symlog', linthresh = 100.0)
				ax.set_ylim( -1.0e+04, 1.0e+04 )
		else:
			ax.plot(X1[:], dset[0,0,:], 'k-', linewidth = 3.25)
			if( UseSymLog ):
				ax.set_yscale('symlog', linthresh = 1.0e0)
				ax.yaxis.set_major_locator(SymmetricalLogLocator(base=10.0, linthresh = 1.0e0, subs=(1.0,)))
				ax.yaxis.set_minor_locator(SymmetricalLogLocator(base=10.0, linthresh = 1.0e0, subs=np.arange(2,10)))
				ax.yaxis.set_major_formatter(LogFormatterMathtext())
				ax.xaxis.set_minor_locator(AutoMinorLocator(5))
			else:
				ax.yaxis.set_minor_locator(AutoMinorLocator(5))
				ax.xaxis.set_minor_locator(AutoMinorLocator(5))
			if (lims != []):
				ax.set_ylim( min(lims), max(lims) )
			else:
				ax.set_ylim( min(dset[0,0,:]), max(dset[0,0,:]) )
			ax.set_title( DirName + "\n" + "Time: " + f"{time[0]:.2f}" + " " + r"$[\mathrm{ms}]$" )
			ax.set_ylabel( Units )
			ax.set_xlabel( r'$\mathrm{r}$' + ' ' + r'$[\mathrm{km}]$' )

	elif dim == 2:

		fig = plt.figure()
		if sys.argv[3] == 'Conserved Magnetic Field (1)':
			plot = plt.pcolormesh(X2, X1, dset[0,:,:], norm = colors.SymLogNorm(linthresh = 1.00e+07, vmin = -1.00e+15, vmax = 1.00e+15, base = 10))
		elif sys.argv[3] == 'Three-Velocity (3)' and sys.argv[3] == 'Sub':
			q = - ( 1.0 / dset[0,:,:] ) * np.gradient( dset[0,:,:], X1, axis = 0, edge_order = 2 )
			for i in range(0, np.shape(X1)[0]):
				q[:,i] = X1[i] * q[:,i]
			gammaMRI = q * dset[0,:,:] / 2.0
			plot = plt.pcolormesh(X2, X1, gammaMRI, vmin = 1.00e+03, vmax = 2.50e+03)
		plt.colorbar()

	iFile = iFile + 1

	print("Writing file # " + str(iFile) + " out of " + str(len(DataFiles)))

	fig.savefig(DirPath + "/{:04}".format(iFile), dpi=200, bbox_inches="tight")

	plt.close(fig)

MovieCommand = "ffmpeg -framerate 30 -i " + os.getcwd() + "/Images_" + DirName + "_" + FieldName + "/%04d.png -vf scale=1024:-2:flags=lanczos,format=yuv420p -r 30 -c:v libx264 -preset slow -crf 20 -tune stillimage -profile:v high -level 4.0 -movflags +faststart -y movie_" + DirName + "_" + FieldName + ".mp4"

os.system( MovieCommand )
