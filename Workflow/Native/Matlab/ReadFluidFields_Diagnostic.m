function...
  [ Time, X1, X2, X3, TCI, Shock_X1, Shock_X2, Shock_X3, Theta_1, Theta_2, Theta_3, Min_E, Max_E ]...
    = ReadFluidFields_Diagnostic( AppName, FileNumber, Directory )

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

  TCI      = h5read( FileName, '/Fluid Fields/Diagnostic/TCI' );
  Shock_X1 = h5read( FileName, '/Fluid Fields/Diagnostic/Shock (X1)' );
  Shock_X2 = h5read( FileName, '/Fluid Fields/Diagnostic/Shock (X2)' );
  Shock_X3 = h5read( FileName, '/Fluid Fields/Diagnostic/Shock (X3)' );
  Theta_1  = h5read( FileName, '/Fluid Fields/Diagnostic/Theta 1' );
  Theta_2  = h5read( FileName, '/Fluid Fields/Diagnostic/Theta 2' );
  Theta_3  = h5read( FileName, '/Fluid Fields/Diagnostic/Theta 3' );
  Min_E    = h5read( FileName, '/Fluid Fields/Diagnostic/Min E' );
  Max_E    = h5read( FileName, '/Fluid Fields/Diagnostic/Max E' );

end
