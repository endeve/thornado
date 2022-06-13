function [ Time, X1, X2, X3, nIter_Outer, nIter_Inner ]...
  = ReadRadiationFields_IterationCounts( AppName, FileNumber, Directory)

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end
  
  FileName = [ DirName '/' AppName '_RadiationFields_' sprintf( '%06d', FileNumber ) '.h5' ];
  
  Time = h5read( FileName, '/Time' );  
  X1   = h5read( FileName, '/Spatial Grid/X1' );
  X2   = h5read( FileName, '/Spatial Grid/X2' );
  X3   = h5read( FileName, '/Spatial Grid/X3' );
  
  nIter_Outer = h5read( FileName, '/Iteration Counts/Outer Iterations' );
  nIter_Inner = h5read( FileName, '/Iteration Counts/Inner Iterations' );

end

