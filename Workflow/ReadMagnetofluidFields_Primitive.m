function...
  [ Time, X1, X2, X3, uPM_D, uPM_V1, uPM_V2, uPM_V3, uPM_E, uPM_Ne,... 
    uPM_B1, uPM_B2, uPM_B3, uPM_Chi]...
    = ReadMagnetofluidFields_Primitive( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_MagnetofluidFields_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  
  X1   = h5read( FileName, '/Spatial Grid/X1' );
  X2   = h5read( FileName, '/Spatial Grid/X2' );
  X3   = h5read( FileName, '/Spatial Grid/X3' );
  
  uPM_D   = h5read( FileName, '/Magnetofluid Fields/Primitive/Comoving Baryon Density' );
  uPM_V1  = h5read( FileName, '/Magnetofluid Fields/Primitive/Three-Velocity (1)' );
  uPM_V2  = h5read( FileName, '/Magnetofluid Fields/Primitive/Three-Velocity (2)' );
  uPM_V3  = h5read( FileName, '/Magnetofluid Fields/Primitive/Three-Velocity (3)' );
  uPM_E   = h5read( FileName, '/Magnetofluid Fields/Primitive/Internal Energy Density' );
  uPM_Ne  = h5read( FileName, '/Magnetofluid Fields/Primitive/Comoving Electron Density' );
  uPM_B1  = h5read( FileName, '/Magnetofluid Fields/Primitive/Eulerian Magnetic Field (1)' );
  uPM_B2  = h5read( FileName, '/Magnetofluid Fields/Primitive/Eulerian Magnetic Field (2)' );
  uPM_B3  = h5read( FileName, '/Magnetofluid Fields/Primitive/Eulerian Magnetic Field (3)' );
  uPM_Chi = h5read( FileName, '/Magnetofluid Fields/Primitive/Divergence Violation Field' ); 
  
end