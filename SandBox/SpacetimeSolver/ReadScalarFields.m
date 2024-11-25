function...
  [ Time, X1, X2, X3, uSF_u, uSF_v ]...
    = ReadScalarFields( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_ScalarWave_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  
  X1   = h5read( FileName, '/Spatial Grid/X1' );
  X2   = h5read( FileName, '/Spatial Grid/X2' );
  X3   = h5read( FileName, '/Spatial Grid/X3' );
  
  uSF_u     = h5read( FileName, '/Scalar Fields/U' );
  uSF_v     = h5read( FileName, '/Scalar Fields/V' );
  
end