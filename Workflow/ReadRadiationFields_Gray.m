function...
  [ Time, X1, X2, X3, uGR_N, uGR_D, uGR_I1, uGR_I2, uGR_I3, uGR_J, uGR_H1, uGR_H2, uGR_H3, uGR_RMS, uGR_F, uGR_K, uGR_Q ]...
    = ReadRadiationFields_Gray( AppName, FileNumber, Species, Directory )

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

  FileName = [ DirName '/' AppName '_RadiationFields_' sprintf( '%06d', FileNumber ) '.h5' ];

  Time = h5read( FileName, '/Time' );
  X1   = h5read( FileName, '/Spatial Grid/X1' );
  X2   = h5read( FileName, '/Spatial Grid/X2' );
  X3   = h5read( FileName, '/Spatial Grid/X3' );
  
  SpeciesString = [ 'Species_' sprintf( '%02d', SpeciesIndex ) ];
  
  uGR_N  = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Eulerian Number Density' ] );

  uGR_D  = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Lagrangian Number Density' ] );
  uGR_I1 = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Lagrangian Number Flux Density (1)' ] );
  uGR_I2 = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Lagrangian Number Flux Density (2)' ] );
  uGR_I3 = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Lagrangian Number Flux Density (3)' ] );
  
  uGR_J  = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Lagrangian Energy Density' ] );
  uGR_H1 = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Lagrangian Energy Flux Density (1)' ] );
  uGR_H2 = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Lagrangian Energy Flux Density (2)' ] );
  uGR_H3 = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Lagrangian Energy Flux Density (3)' ] );

  uGR_RMS = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/RMS Energy' ] );

  uGR_F = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Lagrangian Flux Factor' ] );
  uGR_K = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Lagrangian Eddington Factor' ] );
  uGR_Q = h5read( FileName, [ '/Radiation Fields/' SpeciesString '/Gray/Lagrangian Heat Flux Factor' ] );

end