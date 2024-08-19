function [ T, X1, MinIt, MaxIt, AveIt, AveIt_Inner ]...
  = GetTrace_NonlinearSolverTally( lo, hi, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end
  
  nT = hi - lo + 1;
  
  FileName = [ DirName '/' 'NonlinearSolverTally_' sprintf( '%06d', lo ) '.h5' ];
  
  X1 = h5read( FileName, '/Spatial Grid/X1' );  
  nX = size( X1, 1 );
  
  T = zeros( nT, 1 );
  MinIt = zeros( nT, nX );
  MaxIt = zeros( nT, nX );
  AveIt = zeros( nT, nX );
  AveIt_Inner = zeros( nT, nX );
  
  for i = 1 : nT
    
    [ T(i), ~, ~, ~, MinIt(i,1:nX), MaxIt(i,1:nX), AveIt(i,1:nX), AveIt_Inner(i,1:nX) ]...
      = ReadNonlinearSolverTally( lo+(i-1), DirName );
  
  end

end