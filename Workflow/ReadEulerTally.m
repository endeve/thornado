function [ Time, D, S1, S2, S3, E, E_i, E_k, E_g ]...
  = ReadEulerTally( FileName, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end
  
  FullFileName = [ DirName '/' FileName ];
  
  n = CountLines( FullFileName ) - 1;

  Time = zeros( n, 1 );
  D    = zeros( n, 1 );
  S1   = zeros( n, 1 );
  S2   = zeros( n, 1 );
  S3   = zeros( n, 1 );
  E    = zeros( n, 1 );
  E_i  = zeros( n, 1 );
  E_k  = zeros( n, 1 );
  E_g  = zeros( n, 1 );
  
  fileID = fopen( FullFileName );
  line   = fgetl( fileID );
  offset = 20;
  for i = 1 : n
    line = fgetl( fileID );

    Time(i) = str2double( line(01+00*offset:00+01*offset) );
    D   (i) = str2double( line(01+01*offset:01+02*offset) );
    S1  (i) = str2double( line(02+02*offset:02+03*offset) );
    S2  (i) = str2double( line(03+03*offset:03+04*offset) );
    S3  (i) = str2double( line(04+04*offset:04+05*offset) );
    E   (i) = str2double( line(05+05*offset:05+06*offset) );
    E_i (i) = str2double( line(06+06*offset:06+07*offset) );
    E_k (i) = str2double( line(07+07*offset:07+08*offset) );
    E_g (i) = str2double( line(08+08*offset:08+09*offset) );
    
  end
  fclose( fileID );

end