function...
  [ Time, X1_C, X2_C, X3_C, nIter_Outer, nIter_Inner, PL_Theta_1, PL_Theta_2, PL_dEnergy ]...
    = ReadRadiationFields_Diagnostic( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_RadiationFields_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  X1_C = h5read( FileName, '/Spatial Grid/X1_C' );
  X2_C = h5read( FileName, '/Spatial Grid/X2_C' );
  X3_C = h5read( FileName, '/Spatial Grid/X3_C' );
  
  nIter_Outer = h5read( FileName, '/Radiation Fields/Diagnostic/Outer Iterations' );
  nIter_Inner = h5read( FileName, '/Radiation Fields/Diagnostic/Inner Iterations' );
  PL_Theta_1  = h5read( FileName, '/Radiation Fields/Diagnostic/Positivity Limiter Theta 1' );
  PL_Theta_2  = h5read( FileName, '/Radiation Fields/Diagnostic/Positivity Limiter Theta 2' );
  PL_dEnergy  = h5read( FileName, '/Radiation Fields/Diagnostic/Positivity Limiter Energy Change' );
  
end