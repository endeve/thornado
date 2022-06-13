function [ nIter_Outer, nIter_Inner ]...
  = ReadRadiationFields_IterationCounts( AppName, FileNumber, Directory)

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end
  
  FileName = [ DirName '/' AppName '_RadiationFields_' sprintf( '%06d', FileNumber ) '.h5' ];
  
  nIter_Outer = h5read( FileName, '/Iteration Counts/Outer Iterations' );
  nIter_Inner = h5read( FileName, '/Iteration Counts/Inner Iterations' );

end

