close all
clear all

k = 1;
d = 2;

i = 0;
for iy = 1 : k+1
for ix = 1 : k+1
  i = i + 1;
  indQ(1,i) = ix;
  indQ(2,i) = iy;
end
end

[ xN, ~ ]  = GetQuadrature( 1*(k+1), 'LG' ); % Nodal Points in DG Method;
[ xQ, wQ ] = GetQuadrature( 2*(k+1), 'LG' ); % Quadrature Points and Weights

% --- Compute P (Projection) and R (Reconstruction) Matrices (specialize to k = 1, d = 2) ---

P = zeros( 4, 4 );

for j = 1 : 4
for i = 1 : 4
  
  % --- P ---
  
  switch i
    case 1
      P(i,j) = sum( wQ(:) .* LagrangeP( 0.5.*(xQ(:)-0.5), indQ(1,j), xN, 2 ) )...
             * sum( wQ(:) .* LagrangeP( 0.5.*(xQ(:)-0.5), indQ(2,j), xN, 2 ) );
    case 2
      P(i,j) = sum( wQ(:) .* LagrangeP( 0.5.*(xQ(:)+0.5), indQ(1,j), xN, 2 ) )...
             * sum( wQ(:) .* LagrangeP( 0.5.*(xQ(:)-0.5), indQ(2,j), xN, 2 ) );
    case 3
      P(i,j) = sum( wQ(:) .* LagrangeP( 0.5.*(xQ(:)-0.5), indQ(1,j), xN, 2 ) )...
             * sum( wQ(:) .* LagrangeP( 0.5.*(xQ(:)+0.5), indQ(2,j), xN, 2 ) );
    case 4
      P(i,j) = sum( wQ(:) .* LagrangeP( 0.5.*(xQ(:)+0.5), indQ(1,j), xN, 2 ) )...
             * sum( wQ(:) .* LagrangeP( 0.5.*(xQ(:)+0.5), indQ(2,j), xN, 2 ) );
  end
  
end
end
R = inv( P );

fig_1 = figure( 1 );

clims = [-0.25,1.25];
imagesc( P, clims )
set( gca, 'fontsize',  15 )
set( gca, 'linewidth', 02 )
set( gca, 'TickDir', 'out')
axis square
colorbar
colormap hot
xticks([1 2 3 4])
yticks([1 2 3 4])
title( 'Projection Matrix, P', 'fontsize', 15 )

print( fig_1, '-dpng', './ProjectionMatrix.png' )

fig_2 = figure( 2 );

clims = [-0.25,1.25];
imagesc( R, clims )
set( gca, 'fontsize',  15 )
set( gca, 'linewidth', 02 )
set( gca, 'TickDir', 'out')
axis square
colorbar
colormap hot
xticks([1 2 3 4])
yticks([1 2 3 4])
title( 'Reconstruction Matrix, R', 'fontsize', 15 )

print( fig_2, '-dpng', './ReconstructionMatrix.png' )