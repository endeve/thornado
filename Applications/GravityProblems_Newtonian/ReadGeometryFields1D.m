function...
  [ Time, X, uGF_Phi_N, uGF_Gm_dd_11, uGF_Gm_dd_22, uGF_Gm_dd_33 ]...
    = ReadGeometryFields1D( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_GeometryFields_'...
                 sprintf( '%06d', FileNumber ) '.dat' ];

  GF = load( FileName )';

  iO = 2;

  Time         = GF(1);
  nX           = GF(2);
  X            = GF(iO+1+00*nX:iO+01*nX);
  uGF_Phi_N    = GF(iO+1+01*nX:iO+02*nX);
  uGF_Gm_dd_11 = GF(iO+1+02*nX:iO+03*nX);
  uGF_Gm_dd_22 = GF(iO+1+03*nX:iO+04*nX);
  uGF_Gm_dd_33 = GF(iO+1+04*nX:iO+05*nX);
  
%   keyboard

end

