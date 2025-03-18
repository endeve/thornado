function...
  [ Time, X1, X2, X3, uAM_P, uAM_Pb, uAM_T, uAM_Ye, uAM_S, uAM_E,...
    uAM_h, uAM_hb, uAM_Me, uAM_Mp, uAM_Mn, uAM_Yp, uAM_Yn, uAM_Ya,...
    uAM_Yh, uAM_Gm, uAM_Cs, uAM_Ca, uAM_EF1, uAM_EF2, uAM_EF3,...
    uAM_Tem00, uAM_Tem11, uAM_Tem22, uAM_Tem33,...
    uAM_Tem12, uAM_Tem13, uAM_Tem23 ]...
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
  
  uAM_P     = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Pressure' );
  uAM_Pb    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Magnetic Pressure');
  uAM_T     = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Temperature' );
  uAM_Ye    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Electron Fraction' );
  uAM_S     = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Entropy Per Baryon' );
  uAM_E     = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Specific Internal Energy' );
  uAM_h     = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Specific Enthalpy');
  uAM_hb    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Specific Magnetic Enthalpy');
  uAM_Me    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Electron Chemical Potential' );
  uAM_Mp    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Proton Chemical Potential' );
  uAM_Mn    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Neutron Chemical Potential' );
  uAM_Yp    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Proton Mass Fraction' );
  uAM_Yn    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Neutron Mass Fraction' );
  uAM_Ya    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Alpha Mass Fraction' );
  uAM_Yh    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Heavy Mass Fraction' );
  uAM_Gm    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Ratio of Specific Heats (Gamma)' );
  uAM_Cs    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Sound Speed' );
  uAM_Ca    = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Alfven Speed');
  uAM_EF1   = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Eulerian Electric Field (1)');
  uAM_EF2   = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Eulerian Electric Field (2)');
  uAM_EF3   = h5read( FileName, '/Magnetofluid Fields/Auxiliary/Eulerian Electric Field (3)');
  uAM_Tem00 = h5read( FileName, '/Magnetofluid Fields/Auxiliary/EM Stress Tensor (00)');
  uAM_Tem11 = h5read( FileName, '/Magnetofluid Fields/Auxiliary/EM Stress Tensor (11)');
  uAM_Tem22 = h5read( FileName, '/Magnetofluid Fields/Auxiliary/EM Stress Tensor (22)');
  uAM_Tem33 = h5read( FileName, '/Magnetofluid Fields/Auxiliary/EM Stress Tensor (33)');
  uAM_Tem12 = h5read( FileName, '/Magnetofluid Fields/Auxiliary/EM Stress Tensor (12)');
  uAM_Tem13 = h5read( FileName, '/Magnetofluid Fields/Auxiliary/EM Stress Tensor (13)');
  uAM_Tem23 = h5read( FileName, '/Magnetofluid Fields/Auxiliary/EM Stress Tensor (23)');

end
