function [ Time, Interior, OffGrid, Initial, Change ]...
  = ReadConservationTally( FileName, Directory, C_F )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end
  
  if( ~exist( 'C_F', 'var' ) )
    C_F = 1.0;
  end
  
  FullFileName = [ DirName '/' FileName ];
  
  n = CountLines( FullFileName ) - 1;

  Time     = zeros( n, 1 );
  Interior = zeros( n, 1 );
  OffGrid  = zeros( n, 1 );
  Initial  = zeros( n, 1 );
  Change   = zeros( n, 1 );
  
  fileID = fopen( FullFileName );
  line   = fgetl( fileID );
  offset = 20;
  for i = 1 : n
    line = fgetl( fileID );

    Time    (i) = str2double( line(01+00*offset:00+01*offset) );
    Interior(i) = str2double( line(01+01*offset:01+02*offset) );
    OffGrid (i) = str2double( line(02+02*offset:02+03*offset) );
    Initial (i) = str2double( line(03+03*offset:03+04*offset) );
    Change  (i) = str2double( line(04+04*offset:04+05*offset) );
    
  end
  fclose( fileID );
  
  Interior = C_F * Interior;
  OffGrid  = C_F * OffGrid;
  Initial  = C_F * Initial;
  Change   = C_F * Change;

end