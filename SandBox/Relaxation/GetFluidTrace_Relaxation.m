function [ t, D, T, Ye, Me, Mp, Mn ]...
  = GetFluidTrace_Relaxation( Directory, AppName, lo, hi )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end
  
  nFiles = hi - lo + 1;
  
  t  = zeros( nFiles, 1 );
  D  = zeros( nFiles, 1 );
  T  = zeros( nFiles, 1 );
  Ye = zeros( nFiles, 1 );
  Me = zeros( nFiles, 1 );
  Mp = zeros( nFiles, 1 );
  Mn = zeros( nFiles, 1 );
  
  iFile = 0;
  for i = lo : hi
      
    iFile = iFile + 1;
    
    % --- Primitive ---
    
    [ t(iFile), ~, ~, ~, tmp1 ]...
      = ReadFluidFields_Primitive( AppName, i, DirName );

    D(iFile) = tmp1(1);
    
    % --- Auxiliary ---
    
    [ ~, ~, ~, ~, ~, tmp1, tmp2, ~, ~, tmp3, tmp4, tmp5 ]...
      = ReadFluidFields_Auxiliary( AppName, i, DirName );
  
    T (iFile) = tmp1(1);
    Ye(iFile) = tmp2(1);
    Me(iFile) = tmp3(1);
    Mp(iFile) = tmp4(1);
    Mn(iFile) = tmp5(1);
    
  end

end

