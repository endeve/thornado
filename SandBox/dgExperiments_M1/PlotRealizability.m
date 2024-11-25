function [ ] = PlotRealizability( J, H )

  plot( squeeze( H(:) ), squeeze( J(:) ), '.m' ); hold on
  
  yy_1 = linspace( 0.0, 1.0, 1024 )';
  yy_2 = linspace( 0.0, 2.0, 1024 )';
  xx_1 = (1.0-yy_1).*yy_1;
  xx_2 = yy_2;

  plot( + xx_2, yy_2, '-r', 'linewidth', 2 )
  plot( - xx_2, yy_2, '-r', 'linewidth', 2 )
  plot( + xx_1, yy_1, '-k', 'linewidth', 2 )
  plot( - xx_1, yy_1, '-k', 'linewidth', 2 )
  
  xMin = - 0.6 * max( J(:) );
  xMax = + 0.6 * max( J(:) );
  yMin = - 0.01 * max( [ 1.0 max( J(:) ) ] );
  yMax = + 1.01 * max( [ 1.0 max( J(:) ) ] );
  
  axis( [ xMin xMax yMin yMax ] );
    
  hold off

end

