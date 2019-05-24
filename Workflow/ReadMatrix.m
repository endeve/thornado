function [ Matrix ] = ReadMatrix( FileName )

  TMP = load( FileName )';
  M = TMP(1); N = TMP(2);
  
  Matrix = zeros( M, N );
  for iCOL = 1 : N
    Matrix(:,iCOL)...
      = TMP((iCOL-1)*M+3:iCOL*M+2);
  end

end

