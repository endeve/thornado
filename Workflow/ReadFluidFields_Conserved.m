function...
  [ Time, X1, X2, X3, uCF_D, uCF_S1, uCF_S2, uCF_S3, uCF_E, uCF_Ne ]...
    = ReadFluidFields_Conserved( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_FluidFields_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  
  X1   = h5read( FileName, '/Spatial Grid/X1' );
  X2   = h5read( FileName, '/Spatial Grid/X2' );
  X3   = h5read( FileName, '/Spatial Grid/X3' );
  
  uCF_D  = h5read( FileName, '/Fluid Fields/Conserved/Conserved Baryon Density' );
  uCF_S1 = h5read( FileName, '/Fluid Fields/Conserved/Conserved Momentum Density (1)' );
  uCF_S2 = h5read( FileName, '/Fluid Fields/Conserved/Conserved Momentum Density (2)' );
  uCF_S3 = h5read( FileName, '/Fluid Fields/Conserved/Conserved Momentum Density (3)' );
  uCF_E  = h5read( FileName, '/Fluid Fields/Conserved/Conserved Energy Density' );
  uCF_Ne = h5read( FileName, '/Fluid Fields/Conserved/Conserved Electron Density' );
  
end