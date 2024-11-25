function [ R ] = ComputeRadius( X1, X2, X3 )

  nR = prod( [ size(X1,1), size(X2,1), size(X3,1) ] );
  R  = zeros( nR, 1 );
  
  iR = 0;
  for iX3 = 1 : size( X3, 1 )
  for iX2 = 1 : size( X2, 1 )
  for iX1 = 1 : size( X1, 1 )
    
    iR = iR + 1;
    
    R(iR) = sqrt( X1(iX1)^2 + X2(iX2)^2 + X3(iX3)^2 );
    
  end
  end
  end

end

