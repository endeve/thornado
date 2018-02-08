function...
  [ Time, E, X1, X2, X3, uCR_N, uCR_G1, uCR_G2, uCR_G3, uPR_D, uPR_I1, uPR_I2, uPR_I3 ]...
    = ReadRadiationFields( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_RadiationFields_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  E    = h5read( FileName, '/Energy Grid/E' );
  X1   = h5read( FileName, '/Spatial Grid/X1' );
  X2   = h5read( FileName, '/Spatial Grid/X2' );
  X3   = h5read( FileName, '/Spatial Grid/X3' );
  
  uCR_N  = h5read( FileName, '/Radiation Fields/Species_01/Conserved/Eulerian Number Density' );
  uCR_G1 = h5read( FileName, '/Radiation Fields/Species_01/Conserved/Eulerian Number Flux Density (1)' );
  uCR_G2 = h5read( FileName, '/Radiation Fields/Species_01/Conserved/Eulerian Number Flux Density (2)' );
  uCR_G3 = h5read( FileName, '/Radiation Fields/Species_01/Conserved/Eulerian Number Flux Density (3)' );
  
  uPR_D  = h5read( FileName, '/Radiation Fields/Species_01/Primitive/Lagrangian Number Density' );
  uPR_I1 = h5read( FileName, '/Radiation Fields/Species_01/Primitive/Lagrangian Number Flux Density (1)' );
  uPR_I2 = h5read( FileName, '/Radiation Fields/Species_01/Primitive/Lagrangian Number Flux Density (2)' );
  uPR_I3 = h5read( FileName, '/Radiation Fields/Species_01/Primitive/Lagrangian Number Flux Density (3)' );

end

