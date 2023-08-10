function...
  [ Time, E, X1, X2, X3, uCM_E, uCM_F1, uCM_F2, uCM_F3, uPM_J, uPM_H1, uPM_H2, uPM_H3 ]...
    = ReadTwoMomentFields( AppName, FileNumber, Species, Directory )

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

  FileName = [ DirName '/' AppName '_TwoMomentFields_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  E    = h5read( FileName, '/Energy Grid/E' );
  X1   = h5read( FileName, '/Spatial Grid/X1' );
  X2   = h5read( FileName, '/Spatial Grid/X2' );
  X3   = h5read( FileName, '/Spatial Grid/X3' );
  
  SpeciesString = [ 'Species_' sprintf( '%02d', SpeciesIndex ) ];
  
  uCM_E  = h5read( FileName, [ '/Two-Moment Fields/' SpeciesString '/Conserved/Eulerian Energy Density' ] );
  uCM_F1 = h5read( FileName, [ '/Two-Moment Fields/' SpeciesString '/Conserved/Eulerian Momentum Density (1)' ] );
  uCM_F2 = h5read( FileName, [ '/Two-Moment Fields/' SpeciesString '/Conserved/Eulerian Momentum Density (2)' ] );
  uCM_F3 = h5read( FileName, [ '/Two-Moment Fields/' SpeciesString '/Conserved/Eulerian Momentum Density (3)' ] );
  
  uPM_J  = h5read( FileName, [ '/Two-Moment Fields/' SpeciesString '/Primitive/Lagrangian Energy Density' ] );
  uPM_H1 = h5read( FileName, [ '/Two-Moment Fields/' SpeciesString '/Primitive/Lagrangian Momentum Density (1)' ] );
  uPM_H2 = h5read( FileName, [ '/Two-Moment Fields/' SpeciesString '/Primitive/Lagrangian Momentum Density (2)' ] );
  uPM_H3 = h5read( FileName, [ '/Two-Moment Fields/' SpeciesString '/Primitive/Lagrangian Momentum Density (3)' ] );

end