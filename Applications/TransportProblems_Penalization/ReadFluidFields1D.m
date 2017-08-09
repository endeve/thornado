function...
  [ Time, X, uCF_D, uCF_S1, uCF_S2, uCF_S3, uCF_E, uCF_Ne,...
    uPF_D, uPF_V1, uPF_V2, uPF_V3, uPF_E, uPF_Ne, uAF_P, uAF_T,...
    uAF_Ye, uAF_S, uAF_E, uAF_Me, uAF_Mp, uAF_Mn, uAF_Xp, uAF_Xn,...
    uAF_Xa, uAF_Xh, uAF_Gm, uAF_Cs, Shock ]...
    = ReadFluidFields1D( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_FluidFields_'...
                 sprintf( '%06d', FileNumber ) '.dat' ];

  FF = load( FileName )';

  iO = 2;

  Time   = FF(1);
  nX     = FF(2);
  X      = FF(iO+1+00*nX:iO+01*nX);
  uCF_D  = FF(iO+1+01*nX:iO+02*nX);
  uCF_S1 = FF(iO+1+02*nX:iO+03*nX);
  uCF_S2 = FF(iO+1+03*nX:iO+04*nX);
  uCF_S3 = FF(iO+1+04*nX:iO+05*nX);
  uCF_E  = FF(iO+1+05*nX:iO+06*nX);
  uCF_Ne = FF(iO+1+06*nX:iO+07*nX);
  uPF_D  = FF(iO+1+07*nX:iO+08*nX);
  uPF_V1 = FF(iO+1+08*nX:iO+09*nX);
  uPF_V2 = FF(iO+1+09*nX:iO+10*nX);
  uPF_V3 = FF(iO+1+10*nX:iO+11*nX);
  uPF_E  = FF(iO+1+11*nX:iO+12*nX);
  uPF_Ne = FF(iO+1+12*nX:iO+13*nX);
  uAF_P  = FF(iO+1+13*nX:iO+14*nX);
  uAF_T  = FF(iO+1+14*nX:iO+15*nX);
  uAF_Ye = FF(iO+1+15*nX:iO+16*nX);
  uAF_S  = FF(iO+1+16*nX:iO+17*nX);
  uAF_E  = FF(iO+1+17*nX:iO+18*nX);
  uAF_Me = FF(iO+1+18*nX:iO+19*nX);
  uAF_Mp = FF(iO+1+19*nX:iO+20*nX);
  uAF_Mn = FF(iO+1+20*nX:iO+21*nX);
  uAF_Xp = FF(iO+1+21*nX:iO+22*nX);
  uAF_Xn = FF(iO+1+22*nX:iO+23*nX);
  uAF_Xa = FF(iO+1+23*nX:iO+24*nX);
  uAF_Xh = FF(iO+1+24*nX:iO+25*nX);  
  uAF_Gm = FF(iO+1+25*nX:iO+26*nX);
  uAF_Cs = FF(iO+1+26*nX:iO+27*nX);
  Shock  = FF(iO+1+27*nX:iO+28*nX);
  
  keyboard

end

