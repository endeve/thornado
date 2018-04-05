function...
  [ Time, X1, X2, X3, uAF_P, uAF_T, uAF_Ye, uAF_S, uAF_E, uAF_Me, uAF_Mp, uAF_Mn, uAF_Yp, uAF_Yn, uAF_Ya, uAF_Yh, uAF_Gm, uAF_Cs ]...
    = ReadFluidFields_Auxiliary( AppName, FileNumber, Directory )

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
  
  uAF_P  = h5read( FileName, '/Fluid Fields/Auxiliary/Pressure' );
  uAF_T  = h5read( FileName, '/Fluid Fields/Auxiliary/Temperature' );
  uAF_Ye = h5read( FileName, '/Fluid Fields/Auxiliary/Electron Fraction' );
  uAF_S  = h5read( FileName, '/Fluid Fields/Auxiliary/Entropy Per Baryon' );
  uAF_E  = h5read( FileName, '/Fluid Fields/Auxiliary/Specific Internal Energy' );
  uAF_Me = h5read( FileName, '/Fluid Fields/Auxiliary/Electron Chemical Potential' );
  uAF_Mp = h5read( FileName, '/Fluid Fields/Auxiliary/Proton Chemical Potential' );
  uAF_Mn = h5read( FileName, '/Fluid Fields/Auxiliary/Neutron Chemical Potential' );
  uAF_Yp = h5read( FileName, '/Fluid Fields/Auxiliary/Proton Mass Fraction' );
  uAF_Yn = h5read( FileName, '/Fluid Fields/Auxiliary/Neutron Mass Fraction' );
  uAF_Ya = h5read( FileName, '/Fluid Fields/Auxiliary/Alpha Mass Fraction' );
  uAF_Yh = h5read( FileName, '/Fluid Fields/Auxiliary/Heavy Mass Fraction' );
  uAF_Gm = h5read( FileName, '/Fluid Fields/Auxiliary/Ratio of Specific Heats (Gamma)' );
  uAF_Cs = h5read( FileName, '/Fluid Fields/Auxiliary/Sound Speed' );
  
end