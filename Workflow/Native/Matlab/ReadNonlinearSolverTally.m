function [ t, X1, X2, X3, MinIt, MaxIt, AveIt, AveIt_Inner ]...
  = ReadNonlinearSolverTally( FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' 'NonlinearSolverTally_' sprintf( '%06d', FileNumber ) '.h5' ];
  
  t = h5read( FileName, '/Time' );
  
  X1 = h5read( FileName, '/Spatial Grid/X1' );
  X2 = h5read( FileName, '/Spatial Grid/X2' );
  X3 = h5read( FileName, '/Spatial Grid/X3' );
  
  MinIt       = h5read( FileName, '/Iteration Counts/Min Iterations' );
  MaxIt       = h5read( FileName, '/Iteration Counts/Max Iterations' );
  AveIt       = h5read( FileName, '/Iteration Counts/Average Iterations' );
  AveIt_Inner = h5read( FileName, '/Iteration Counts/Average Iterations (Inner)' );
  
end