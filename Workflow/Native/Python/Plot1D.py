import h5py as h5
import matplotlib.pyplot as plt

def PlotField1D():

  File = '../../../SandBox/dgExperiments_MHD_Relativistic_IDEAL/Output/Advection1D_MagnetofluidFields_000101.h5'

  fH5 = h5.File( File, 'r' );

  dsetX  = fH5[ 'Spatial Grid'                  ]
  dsetFF = fH5[ 'Magnetofluid Fields/Conserved' ]

  X1 = dsetX['X1'][()]
  X2 = dsetX['X2'][()]
  X3 = dsetX['X3'][()]

  uCF_D  = dsetFF[ 'Conserved Baryon Density'       ][()]
  uCF_S1 = dsetFF[ 'Conserved Momentum Density (1)' ][()]
  uCF_S2 = dsetFF[ 'Conserved Momentum Density (2)' ][()]
  uCF_S3 = dsetFF[ 'Conserved Momentum Density (3)' ][()]
  uCF_E  = dsetFF[ 'Conserved Energy Density'       ][()]
  uCF_Ne = dsetFF[ 'Conserved Electron Density'     ][()]

  plt.plot(X1, uCF_D[0,0,:])

  plt.show()

  fH5.close()

PlotField1D()
