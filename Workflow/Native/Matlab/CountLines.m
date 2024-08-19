function [ nLines ] = CountLines( FileName )

  fid = fopen( FileName );
  nLines = 0;
  line   = fgetl( fid );
  while ischar( line )
    nLines = nLines + 1;
    line   = fgetl( fid );
  end
  fclose(fid);

end

