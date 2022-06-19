function [ W ] = CreateMesh( L, R, N, Z )
  
  W = zeros( N, 1 );

  if( Z == 1.0 )
    
    W(:) = ( R - L ) / N;
      
  else
      
    W(1) = ( R - L ) * ( Z - 1.0 ) / ( Z^N - 1.0 );
    for i = 2 : N
      
      W(i) = Z * W(i-1);
      
    end
      
  end

end

