function...
  [ Time, X1, X2, X3, uGF_H1, uGF_H2, uGF_H3, uGF_Gm11, uGF_Gm22, uGF_Gm33, uGF_SqrtGm, uGF_Phi_N ]...
    = ReadGeometryFields( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_GeometryFields_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  
  X1   = h5read( FileName, '/Spatial Grid/X1' );
  X2   = h5read( FileName, '/Spatial Grid/X2' );
  X3   = h5read( FileName, '/Spatial Grid/X3' );
  
  uGF_H1     = h5read( FileName, '/Geometry Fields/Spatial Scale Factor (1)' );
  uGF_H2     = h5read( FileName, '/Geometry Fields/Spatial Scale Factor (2)' );
  uGF_H3     = h5read( FileName, '/Geometry Fields/Spatial Scale Factor (3)' );

  uGF_Gm11   = h5read( FileName, '/Geometry Fields/Spatial Metric Component (11)' );
  uGF_Gm22   = h5read( FileName, '/Geometry Fields/Spatial Metric Component (22)' );
  uGF_Gm33   = h5read( FileName, '/Geometry Fields/Spatial Metric Component (33)' );
  uGF_SqrtGm = h5read( FileName, '/Geometry Fields/Sqrt of Spatial Metric Determina' );

  uGF_Phi_N  = h5read( FileName, '/Geometry Fields/Newtonian Potential' );
  
end