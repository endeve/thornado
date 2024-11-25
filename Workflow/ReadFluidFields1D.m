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

  iO = 3;

  Time   = FF(1);
  nX     = FF(2);
  nNX    = FF(3);
  
  nX_G   = nX * nNX;
  
  X      = FF(iO+1+00*nX_G:iO+01*nX_G);
  uCF_D  = FF(iO+1+01*nX_G:iO+02*nX_G);
  uCF_S1 = FF(iO+1+02*nX_G:iO+03*nX_G);
  uCF_S2 = FF(iO+1+03*nX_G:iO+04*nX_G);
  uCF_S3 = FF(iO+1+04*nX_G:iO+05*nX_G);
  uCF_E  = FF(iO+1+05*nX_G:iO+06*nX_G);
  uCF_Ne = FF(iO+1+06*nX_G:iO+07*nX_G);
  uPF_D  = FF(iO+1+07*nX_G:iO+08*nX_G);
  uPF_V1 = FF(iO+1+08*nX_G:iO+09*nX_G);
  uPF_V2 = FF(iO+1+09*nX_G:iO+10*nX_G);
  uPF_V3 = FF(iO+1+10*nX_G:iO+11*nX_G);
  uPF_E  = FF(iO+1+11*nX_G:iO+12*nX_G);
  uPF_Ne = FF(iO+1+12*nX_G:iO+13*nX_G);
  uAF_P  = FF(iO+1+13*nX_G:iO+14*nX_G);
  uAF_T  = FF(iO+1+14*nX_G:iO+15*nX_G);
  uAF_Ye = FF(iO+1+15*nX_G:iO+16*nX_G);
  uAF_S  = FF(iO+1+16*nX_G:iO+17*nX_G);
  uAF_E  = FF(iO+1+17*nX_G:iO+18*nX_G);
  uAF_Me = FF(iO+1+18*nX_G:iO+19*nX_G);
  uAF_Mp = FF(iO+1+19*nX_G:iO+20*nX_G);
  uAF_Mn = FF(iO+1+20*nX_G:iO+21*nX_G);
  uAF_Xp = FF(iO+1+21*nX_G:iO+22*nX_G);
  uAF_Xn = FF(iO+1+22*nX_G:iO+23*nX_G);
  uAF_Xa = FF(iO+1+23*nX_G:iO+24*nX_G);
  uAF_Xh = FF(iO+1+24*nX_G:iO+25*nX_G);  
  uAF_Gm = FF(iO+1+25*nX_G:iO+26*nX_G);
  uAF_Cs = FF(iO+1+26*nX_G:iO+27*nX_G);
  
  Shock = zeros(nNX,1);
  for i = 1 : nX
    Shock((i-1)*nNX+1:i*nNX)...
      = FF(iO+27*nX_G+i:iO+27*nX_G+i);
  end
  
  keyboard

end

