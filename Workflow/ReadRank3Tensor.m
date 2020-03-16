function [ Rank3Tensor ] = ReadRank3Tensor( FileName )

  TMP = load( FileName )';
  M = TMP(1); N = TMP(2); P = TMP(3);
    
  Rank3Tensor = zeros( M, N, P );
  for i3 = 1 : P
    for i2 = 1 : N
      Rank3Tensor(:,i2,i3)...
        = TMP( ( 4 + (i2 - 1)*M + (i3 - 1)*M*N ) : ...
               ( 3 + i2*M + (i3 - 1)*M*N ) );
    end
  end
    
end
        
