function...
  [ Time, X1, X2, X3, uCM_D, uCM_S1, uCM_S2, uCM_S3, uCM_E, uCM_Ne,...
    uCM_B1, uCM_B2, uCM_B3, uCM_Chi]...
    = ReadMagnetofluidFields_Conserved( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_MagnetofluidFields_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  
  X1   = h5read( FileName, '/Spatial Grid/X1' );
  X2   = h5read( FileName, '/Spatial Grid/X2' );
  X3   = h5read( FileName, '/Spatial Grid/X3' );
  
  uCM_D   = h5read( FileName, '/Magnetofluid Fields/Conserved/Conserved Baryon Density' );
  uCM_S1  = h5read( FileName, '/Magnetofluid Fields/Conserved/Conserved Momentum Density (1)' );
  uCM_S2  = h5read( FileName, '/Magnetofluid Fields/Conserved/Conserved Momentum Density (2)' );
  uCM_S3  = h5read( FileName, '/Magnetofluid Fields/Conserved/Conserved Momentum Density (3)' );
  uCM_E   = h5read( FileName, '/Magnetofluid Fields/Conserved/Conserved Energy Density' );
  uCM_Ne  = h5read( FileName, '/Magnetofluid Fields/Conserved/Conserved Electron Density' );
  uCM_B1  = h5read( FileName, '/Magnetofluid Fields/Conserved/Conserved Magnetic Field (1)' );
  uCM_B2  = h5read( FileName, '/Magnetofluid Fields/Conserved/Conserved Magnetic Field (2)' );
  uCM_B3  = h5read( FileName, '/Magnetofluid Fields/Conserved/Conserved Magnetic Field (3)' );
  uCM_Chi = h5read( FileName, '/Magnetofluid Fields/Conserved/Divergence Violation Field' );
end