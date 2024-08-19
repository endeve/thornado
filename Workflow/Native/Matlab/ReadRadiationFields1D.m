function...
  [ Time, E, X, uCR_N, uCR_G1, uPR_D, uPR_I1, Shock ]...
    = ReadRadiationFields1D( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_RadiationFields_'...
                 sprintf( '%06d', FileNumber ) '.dat' ];

  RF = load( FileName )';

  Time = RF(1);
  nE   = RF(2);
  nNE  = RF(3);
  E    = RF(4:nE*nNE+3);
  nX   = RF(nE*nNE+4:nE*nNE+4);
  nNX  = RF(nE*nNE+5:nE*nNE+5);
  X    = RF(nE*nNE+6:nE*nNE+6+nX*nNX-1);

  nE_G = nE * nNE;
  nX_G = nX * nNX;
  
  iO = nE_G + nX_G + 5;
  
  uCR_N  = zeros( nE_G, nX_G );
  uCR_G1 = zeros( nE_G, nX_G );
  uCR_G2 = zeros( nE_G, nX_G );
  uCR_G3 = zeros( nE_G, nX_G );
  uPR_D  = zeros( nE_G, nX_G );
  uPR_I1 = zeros( nE_G, nX_G );
  uPR_I2 = zeros( nE_G, nX_G );
  uPR_I3 = zeros( nE_G, nX_G );
  
  for i = 1 : nX_G
  
    uCR_N (:,i) = RF(iO+0*nE_G*nX_G+(i-1)*nE_G+1:iO+0*nE_G*nX_G+i*nE_G);
    uCR_G1(:,i) = RF(iO+1*nE_G*nX_G+(i-1)*nE_G+1:iO+1*nE_G*nX_G+i*nE_G);
    uCR_G2(:,i) = RF(iO+2*nE_G*nX_G+(i-1)*nE_G+1:iO+2*nE_G*nX_G+i*nE_G);
    uCR_G3(:,i) = RF(iO+3*nE_G*nX_G+(i-1)*nE_G+1:iO+3*nE_G*nX_G+i*nE_G);
    uPR_D (:,i) = RF(iO+4*nE_G*nX_G+(i-1)*nE_G+1:iO+4*nE_G*nX_G+i*nE_G);
    uPR_I1(:,i) = RF(iO+5*nE_G*nX_G+(i-1)*nE_G+1:iO+5*nE_G*nX_G+i*nE_G);
    uPR_I2(:,i) = RF(iO+6*nE_G*nX_G+(i-1)*nE_G+1:iO+6*nE_G*nX_G+i*nE_G);
    uPR_I3(:,i) = RF(iO+7*nE_G*nX_G+(i-1)*nE_G+1:iO+7*nE_G*nX_G+i*nE_G);

  end
  
  iO = iO+8*nE_G*nX_G;
  
  Shock  = zeros( nE_G, nX_G );
  for j = 1 : nX
    for i = 1 : nE
      Shock((i-1)*nNE+1:i*nNE,(j-1)*nNX+1:j*nNX)...
        = RF(iO+(j-1)*nE+i:iO+(j-1)*nE+i);
    end
  end
  
%   keyboard

end

