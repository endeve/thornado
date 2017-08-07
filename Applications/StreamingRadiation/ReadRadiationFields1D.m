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
  E    = RF(3:nE+2);
  nX   = RF(nE+3:nE+3);
  X    = RF(nE+4:nE+4+nX-1);
  
  iO = nE + nX + 3;
  
  uCR_N  = zeros( nE, nX );
  uCR_G1 = zeros( nE, nX );
  uCR_G2 = zeros( nE, nX );
  uCR_G3 = zeros( nE, nX );
  uPR_D  = zeros( nE, nX );
  uPR_I1 = zeros( nE, nX );
  uPR_I2 = zeros( nE, nX );
  uPR_I3 = zeros( nE, nX );
  Shock  = zeros( nE, nX );
  
  for i = 1 : nX
      
    uCR_N (:,i) = RF(iO+0*nE*nX+(i-1)*nE+1:iO+0*nE*nX+i*nE);
    uCR_G1(:,i) = RF(iO+1*nE*nX+(i-1)*nE+1:iO+1*nE*nX+i*nE);
    uCR_G2(:,i) = RF(iO+2*nE*nX+(i-1)*nE+1:iO+2*nE*nX+i*nE);
    uCR_G3(:,i) = RF(iO+3*nE*nX+(i-1)*nE+1:iO+3*nE*nX+i*nE);
    uPR_D (:,i) = RF(iO+4*nE*nX+(i-1)*nE+1:iO+4*nE*nX+i*nE);
    uPR_I1(:,i) = RF(iO+5*nE*nX+(i-1)*nE+1:iO+5*nE*nX+i*nE);
    uPR_I2(:,i) = RF(iO+6*nE*nX+(i-1)*nE+1:iO+6*nE*nX+i*nE);
    uPR_I3(:,i) = RF(iO+7*nE*nX+(i-1)*nE+1:iO+7*nE*nX+i*nE);
    Shock (:,i) = RF(iO+8*nE*nX+(i-1)*nE+1:iO+8*nE*nX+i*nE);

  end
  
  keyboard

end

