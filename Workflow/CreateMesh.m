function [ W, C, E ] = CreateMesh( L, R, N, Z )
  
  W = zeros( N  , 1 );
  C = zeros( N  , 1 );
  E = zeros( N+1, 1 );

  if( Z == 1.0 )
    
    W(:) = ( R - L ) / N;
      
  else
      
    W(1) = ( R - L ) * ( Z - 1.0 ) / ( Z^N - 1.0 );
    for i = 2 : N
      
      W(i) = Z * W(i-1);
      
    end
      
  end

  E(1) = L;
  for i = 2 : N+1;
    E(i) = E(i-1) + W(i);
  end

  for i = 1 : N
    C(i) = 0.5 * ( E(i) + E(i+1) );
  end

end

