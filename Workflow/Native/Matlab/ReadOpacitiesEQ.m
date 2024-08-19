function...
  [ Time, E, X1, X2, X3, f_EQ ]...
    = ReadOpacitiesEQ( AppName, FileNumber, Species, Directory )

  if( exist( 'Species', 'var' ) )
    SpeciesIndex = Species;
  else
    SpeciesIndex = 1;
  end
  
  if( exist( 'Directory', 'var' ) )
    DirName = Directory;
  else
    DirName = './Output';
  end

  FileName = [ DirName '/' AppName '_Opacities_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  E    = h5read( FileName, '/Energy Grid/E' );
  X1   = h5read( FileName, '/Spatial Grid/X1' );
  X2   = h5read( FileName, '/Spatial Grid/X2' );
  X3   = h5read( FileName, '/Spatial Grid/X3' );
  
  SpeciesString = [ 'Species_' sprintf( '%02d', SpeciesIndex ) ];
  
  f_EQ = h5read( FileName, [ '/Opacities/' SpeciesString '/Equilibrium Distribution' ] );

end