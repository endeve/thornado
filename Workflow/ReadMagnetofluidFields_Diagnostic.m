function...
  [ Time, X1, X2, X3, TCI, Shock_X1, Shock_X2, Shock_X3, Theta_1, Theta_2, Theta_3, Min_E, Max_E ]...
    = ReadMagnetofluidFields_Diagnostic( AppName, FileNumber, Directory )

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

  TCI      = h5read( FileName, '/Magnetofluid Fields/Diagnostic/TCI' );
  Shock_X1 = h5read( FileName, '/Magnetofluid Fields/Diagnostic/Shock (X1)' );
  Shock_X2 = h5read( FileName, '/Magnetofluid Fields/Diagnostic/Shock (X2)' );
  Shock_X3 = h5read( FileName, '/Magnetofluid Fields/Diagnostic/Shock (X3)' );
  Theta_1  = h5read( FileName, '/Magnetofluid Fields/Diagnostic/Theta 1' );
  Theta_2  = h5read( FileName, '/Magnetofluid Fields/Diagnostic/Theta 2' );
  Theta_3  = h5read( FileName, '/Magnetofluid Fields/Diagnostic/Theta 3' );
  Min_E    = h5read( FileName, '/Magnetofluid Fields/Diagnostic/Min E' );
  Max_E    = h5read( FileName, '/Magnetofluid Fields/Diagnostic/Max E' );

end
