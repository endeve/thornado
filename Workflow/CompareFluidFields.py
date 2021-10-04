
import h5py  as h5
import numpy as np

def CompareFluidFields( File1, File2 ):
    
    print( "--------------------------------------" )
    print( "--- Compare Fluid Fields -------------" )
    print( "--------------------------------------" )
    
    # --- Read Data From File1 ---
    
    fH5 = h5.File( File1, 'r' );
    
    dsetFF = fH5[ 'Fluid Fields/Conserved' ];
    
    uCF1_D  = dsetFF[ 'Conserved Baryon Density'       ][()];
    uCF1_S1 = dsetFF[ 'Conserved Momentum Density (1)' ][()];
    uCF1_S2 = dsetFF[ 'Conserved Momentum Density (2)' ][()];
    uCF1_S3 = dsetFF[ 'Conserved Momentum Density (3)' ][()];
    uCF1_E  = dsetFF[ 'Conserved Energy Density'       ][()];
    uCF1_Ne = dsetFF[ 'Conserved Electron Density'     ][()];
    
    fH5.close()
    
    # --- Read Data From File2 ---
    
    fH5 = h5.File( File2, 'r' );
    
    dsetFF = fH5[ 'Fluid Fields/Conserved' ];
    
    uCF2_D  = dsetFF[ 'Conserved Baryon Density'       ][()];
    uCF2_S1 = dsetFF[ 'Conserved Momentum Density (1)' ][()];
    uCF2_S2 = dsetFF[ 'Conserved Momentum Density (2)' ][()];
    uCF2_S3 = dsetFF[ 'Conserved Momentum Density (3)' ][()];
    uCF2_E  = dsetFF[ 'Conserved Energy Density'       ][()];
    uCF2_Ne = dsetFF[ 'Conserved Electron Density'     ][()];
    
    fH5.close()
    
    print( '  max|dD|  = ', np.amax( np.abs( uCF2_D .flatten() - uCF1_D .flatten() ) ) )
    print( '  max|dS1| = ', np.amax( np.abs( uCF2_S1.flatten() - uCF1_S1.flatten() ) ) )
    print( '  max|dS2| = ', np.amax( np.abs( uCF2_S2.flatten() - uCF1_S2.flatten() ) ) )
    print( '  max|dS3| = ', np.amax( np.abs( uCF2_S3.flatten() - uCF1_S3.flatten() ) ) )
    print( '  max|dE|  = ', np.amax( np.abs( uCF2_E .flatten() - uCF1_E .flatten() ) ) )
    print( '  max|dNe| = ', np.amax( np.abs( uCF2_Ne.flatten() - uCF1_Ne.flatten() ) ) )
    
    print( "--------------------------------------" )