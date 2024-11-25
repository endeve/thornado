function...
  [ Time, X1, X2, X3, uPF_D, uPF_V1, uPF_V2, uPF_V3, uPF_E, uPF_Ne ]...
    = ReadFluidFields_Primitive( AppName, FileNumber, Directory )

  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_FluidFields_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  
  X1   = h5read( FileName, '/Spatial Grid/X1' );
  X2   = h5read( FileName, '/Spatial Grid/X2' );
  X3   = h5read( FileName, '/Spatial Grid/X3' );
  
  uPF_D  = h5read( FileName, '/Fluid Fields/Primitive/Comoving Baryon Density' );
  uPF_V1 = h5read( FileName, '/Fluid Fields/Primitive/Three-Velocity (1)' );
  uPF_V2 = h5read( FileName, '/Fluid Fields/Primitive/Three-Velocity (2)' );
  uPF_V3 = h5read( FileName, '/Fluid Fields/Primitive/Three-Velocity (3)' );
  uPF_E  = h5read( FileName, '/Fluid Fields/Primitive/Internal Energy Density' );
  uPF_Ne = h5read( FileName, '/Fluid Fields/Primitive/Comoving Electron Density' );
  
end