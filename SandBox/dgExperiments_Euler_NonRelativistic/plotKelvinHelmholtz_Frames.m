close all
clear all

OutputDir = 'Output_KH_128Sq_noIndicator';
nFrames = 100;

for iFrame = 0 : nFrames

  fprintf( '  frame = %d\n', iFrame );

  close all

  % Density: 
  
  [ t, x1, x2, ~, D ]...
    = ReadFluidFields_Conserved...
        ( 'KelvinHelmholtz', iFrame, OutputDir );
  
  fig_1 = figure( 1 );
  set( gcf, 'Visible', 'off' );

  clims = [ 1.0 2.0 ];
  imagesc( x1, x2, D', clims )
  axis( [ 0 1 0 1 ] )
  colormap bone
  title( [ 'Density. t = ' sprintf( '%05d', t ) ] )
  xlabel( 'X', 'fontsize', 20 )
  ylabel( 'Y', 'fontsize', 20 )
  set( gca, 'Ydir', 'normal' )
  colorbar( 'eastoutside' )
  
  print( fig_1, '-dpng',...
         [ OutputDir '/Figures/KelvinHelmholtz_Density_' sprintf( '%05d', iFrame ) '.png' ] )
  
  % Troubled Cell Indicator:
  
  [ t, x1_c, x2_c, ~, Shock ]...
    = ReadShockDetector_Fluid...
        ( 'KelvinHelmholtz', iFrame, OutputDir );
  
  fig_2 = figure( 2 );
  set( gcf, 'Visible', 'off' );

  clims = [ 0.0 0.3 ];
  imagesc( x1_c, x2_c, Shock', clims )
  axis( [ 0 1 0 1 ] )
  colormap bone
  title( [ 'Troubled Cell Indicator. t = ' sprintf( '%05d', t ) ] )
  xlabel( 'X', 'fontsize', 20 )
  ylabel( 'Y', 'fontsize', 20 )
  set( gca, 'Ydir', 'normal' )
  colorbar( 'eastoutside' )
  
  print( fig_2, '-dpng',...
         [ OutputDir '/Figures/KelvinHelmholtz_TroubledCell_' sprintf( '%05d', iFrame ) '.png' ] )
  
end