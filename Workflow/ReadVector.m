function [ Vector ] = ReadVector( FileName )

  TMP = load( FileName )';
  N = TMP(1);
  
  Vector = zeros( N, 1 );

  Vector(1:N) = TMP(2:N+1);

end

