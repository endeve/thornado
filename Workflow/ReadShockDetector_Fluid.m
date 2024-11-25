function...
  [ Time, X1_C, X2_C, X3_C, Shock ]...
    = ReadShockDetector_Fluid( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_FluidFields_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  
  X1_C  = h5read( FileName, '/Spatial Grid/X1_C' );
  X2_C  = h5read( FileName, '/Spatial Grid/X2_C' );
  X3_C  = h5read( FileName, '/Spatial Grid/X3_C' );
  
  Shock = h5read( FileName, '/Shock Detector/Shock' );
  
end