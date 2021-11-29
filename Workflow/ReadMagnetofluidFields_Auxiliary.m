function...
  [ Time, X1, X2, X3, uAM_P, uAM_T, uAM_Ye, uAM_S, uAM_E, uAM_Me, uAM_Mp, uAM_Mn, uAM_Yp, uAM_Yn, uAM_Ya, uAM_Yh, uAM_Gm, uAM_Cs ]...
    = ReadMagnetofluidFields_Auxiliary( AppName, FileNumber, Directory )

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
  
  uAM_P  = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Pressure' );
  uAM_T  = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Temperature' );
  uAM_Ye = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Electron Fraction' );
  uAM_S  = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Entropy Per Baryon' );
  uAM_E  = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Specific Internal Energy' );
  uAM_Me = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Electron Chemical Potential' );
  uAM_Mp = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Proton Chemical Potential' );
  uAM_Mn = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Neutron Chemical Potential' );
  uAM_Yp = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Proton Mass Fraction' );
  uAM_Yn = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Neutron Mass Fraction' );
  uAM_Ya = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Alpha Mass Fraction' );
  uAM_Yh = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Heavy Mass Fraction' );
  uAM_Gm = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Ratio of Specific Heats (Gamma)' );
  uAM_Cs = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Sound Speed' );
  
end