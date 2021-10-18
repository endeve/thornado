
import h5py  as h5
import numpy as np

def CompareRadiationFields( File1, File2 ):
    
    print( "--------------------------------------" )
    print( "--- Compare Radiation Fields ---------" )
    print( "--------------------------------------" )
    
    # --- Read Data From File1 ---
    
    fH5 = h5.File( File1, 'r' );
    
    dsetRF = fH5[ 'Radiation Fields/Species_01/Conserved' ];
    
    uCR1_N  = dsetRF[ 'Eulerian Number Density'          ][()];
    uCR1_G1 = dsetRF[ 'Eulerian Number Flux Density (1)' ][()];
    uCR1_G2 = dsetRF[ 'Eulerian Number Flux Density (2)' ][()];
    uCR1_G3 = dsetRF[ 'Eulerian Number Flux Density (3)' ][()];
    
    fH5.close()
    
    # --- Read Data From File2 ---
    
    fH5 = h5.File( File2, 'r' );
    
    dsetRF = fH5[ 'Radiation Fields/Species_01/Conserved' ];
    
    uCR2_N  = dsetRF[ 'Eulerian Number Density'          ][()];
    uCR2_G1 = dsetRF[ 'Eulerian Number Flux Density (1)' ][()];
    uCR2_G2 = dsetRF[ 'Eulerian Number Flux Density (2)' ][()];
    uCR2_G3 = dsetRF[ 'Eulerian Number Flux Density (3)' ][()];
    
    fH5.close()
    
    print( '  max|dN|  = ', np.amax( np.abs( uCR2_N .flatten() - uCR1_N .flatten() ) ) )
    print( '  max|dG1| = ', np.amax( np.abs( uCR2_G1.flatten() - uCR1_G1.flatten() ) ) )
    print( '  max|dG2| = ', np.amax( np.abs( uCR2_G2.flatten() - uCR1_G2.flatten() ) ) )
    print( '  max|dG3| = ', np.amax( np.abs( uCR2_G3.flatten() - uCR1_G3.flatten() ) ) )
    
    print( "--------------------------------------" )