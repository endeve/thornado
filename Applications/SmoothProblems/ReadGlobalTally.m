function [ Time, MaxD, M, dM, E_F, E_Fi, E_Fk, dE_F,...
           N_F, dN_F, E_G, dE_G, N_R, dN_R, E_R, dE_R ]...
  = ReadGlobalTally( FileName, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end
  
  FullFileName = [ DirName '/' FileName ];
  
  n = CountLines( FullFileName ) - 1;

  Time = zeros( n, 1 );
  MaxD = zeros( n, 1 );
  M    = zeros( n, 1 );
  dM   = zeros( n, 1 );
  E_F  = zeros( n, 1 );
  E_Fi = zeros( n, 1 );
  E_Fk = zeros( n, 1 );
  dE_F = zeros( n, 1 );
  N_F  = zeros( n, 1 );
  dN_F = zeros( n, 1 );
  E_G  = zeros( n, 1 );
  dE_G = zeros( n, 1 );
  N_R  = zeros( n, 1 );
  dN_R = zeros( n, 1 );
  E_R  = zeros( n, 1 );
  dE_R = zeros( n, 1 );
  
  fileID = fopen( FullFileName );
  line   = fgetl( fileID );
  offset = 14;
  for i = 1 : n
    line = fgetl( fileID );

    Time(i) = str2double( line(01+00*offset:00+01*offset) );
    MaxD(i) = str2double( line(01+01*offset:01+02*offset) );
    M   (i) = str2double( line(02+02*offset:02+03*offset) );
    dM  (i) = str2double( line(03+03*offset:03+04*offset) );
    E_F (i) = str2double( line(04+04*offset:04+05*offset) );
    E_Fi(i) = str2double( line(05+05*offset:05+06*offset) );
    E_Fk(i) = str2double( line(06+06*offset:06+07*offset) );
    dE_F(i) = str2double( line(07+07*offset:07+08*offset) );
    N_F (i) = str2double( line(08+08*offset:08+09*offset) );
    dN_F(i) = str2double( line(09+09*offset:09+10*offset) );
    E_G (i) = str2double( line(10+10*offset:10+11*offset) );
    dE_G(i) = str2double( line(11+11*offset:11+12*offset) );
    N_R (i) = str2double( line(12+12*offset:12+13*offset) );
    dN_R(i) = str2double( line(13+13*offset:13+14*offset) );
    E_R (i) = str2double( line(14+14*offset:14+15*offset) );
    dE_R(i) = str2double( line(15+15*offset:15+16*offset) );
    
  end
  fclose( fileID );

end