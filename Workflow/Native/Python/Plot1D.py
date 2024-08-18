import h5py as h5
import matplotlib.pyplot as plt

def PlotField1D():

  File_I = '../../../SandBox/dgExperiments_MHD_Relativistic_IDEAL/Output/Advection1D_MagnetofluidFields_000000.h5'
  File_F = '../../../SandBox/dgExperiments_MHD_Relativistic_IDEAL/Output/Advection1D_MagnetofluidFields_000101.h5'

  fH5_I = h5.File( File_I, 'r' )
  fH5_F = h5.File( File_F, 'r' )

  dsetX_I  = fH5_I   [ 'Spatial Grid'                  ]
  dsetFF_I = fH5_I   [ 'Magnetofluid Fields/Conserved' ]
  X1_I     = dsetX_I [ 'X1'                       ][()]
  uCF_D_I  = dsetFF_I[ 'Conserved Baryon Density' ][()]

  dsetX_F  = fH5_F   [ 'Spatial Grid'                  ]
  dsetFF_F = fH5_F   [ 'Magnetofluid Fields/Conserved' ]
  X1_F     = dsetX_F [ 'X1'                       ][()]
  uCF_D_F  = dsetFF_F[ 'Conserved Baryon Density' ][()]

  plt.plot(X1_I, uCF_D_I[0,0,:])
  plt.plot(X1_F, uCF_D_F[0,0,:])

  plt.legend( [ 'Initial', 'Final' ] )

  plt.xlabel( 'X1'  )
  plt.ylabel( 'rho' )

  plt.show()

  fH5_I.close()
  fH5_F.close()

PlotField1D()
