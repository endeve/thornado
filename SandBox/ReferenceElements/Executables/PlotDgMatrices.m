close all
clear all

% --- Interpolation from Lobatto to Gaussian Points ---

fig_1 = figure( 1 );

L_L2G = ReadMatrix('L_L2G.dat');
clims = [ - 4 0 ];
imagesc( log10( abs( L_L2G ) ), clims )
axis square
title( 'Interpolation Matrix, log_{10}( | L_{ik}^{L2G} | )', 'fontsize', 15 )
xlabel( 'N = 64', 'fontsize', 20 )
ylabel( 'N = 64', 'fontsize', 20 )
colorbar( 'eastoutside' )
colormap gray
print( fig_1, '-dpng', 'L_L2G.png' );

% --- Interpolation to Element Interfaces ---

fig_2 = figure( 2 );

subplot( 3, 2, 1 )

L1L_4 = ReadMatrix('L1L_4.dat');
clims = [ - 1.5 1.5 ];
imagesc( L1L_4, clims )
axis equal
axis( [ 0.5 64.5 0.5 16.5 ] )
title( 'Interpolation Matrix, L_{kl}^{1}(x_{L}^{1})', 'fontsize', 7 )
xlabel( 'N = 64', 'fontsize', 8 )
ylabel( 'N^{1} = 16', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap gray

subplot( 3, 2, 2 )

L1H_4 = ReadMatrix('L1H_4.dat');
clims = [ - 1.5 1.5 ];
imagesc( L1H_4, clims )
axis equal
axis( [ 0.5 64.5 0.5 16.5 ] )
title( 'Interpolation Matrix, L_{kl}^{1}(x_{H}^{1})', 'fontsize', 7 )
xlabel( 'N = 64', 'fontsize', 8 )
ylabel( 'N^{1} = 16', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

subplot( 3, 2, 3 )

L2L_4 = ReadMatrix('L2L_4.dat');
clims = [ - 1.5 1.5 ];
imagesc( L2L_4, clims )
axis equal
axis( [ 0.5 64.5 0.5 16.5 ] )
title( 'Interpolation Matrix, L_{kl}^{2}(x_{L}^{2})', 'fontsize', 7 )
xlabel( 'N = 64', 'fontsize', 8 )
ylabel( 'N^{2} = 16', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

subplot( 3, 2, 4 )

L2H_4 = ReadMatrix('L2H_4.dat');
clims = [ - 1.5 1.5 ];
imagesc( L2H_4, clims )
axis equal
axis( [ 0.5 64.5 0.5 16.5 ] )
title( 'Interpolation Matrix, L_{kl}^{2}(x_{H}^{2})', 'fontsize', 7 )
xlabel( 'N = 64', 'fontsize', 8 )
ylabel( 'N^{2} = 16', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

subplot( 3, 2, 5 )

L3L_4 = ReadMatrix('L3L_4.dat');
clims = [ - 1.5 1.5 ];
imagesc( L3L_4, clims )
axis equal
axis( [ 0.5 64.5 0.5 16.5 ] )
title( 'Interpolation Matrix, L_{kl}^{3}(x_{L}^{3})', 'fontsize', 7 )
xlabel( 'N = 64', 'fontsize', 8 )
ylabel( 'N^{3} = 16', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

subplot( 3, 2, 6 )

L3H_4 = ReadMatrix('L3H_4.dat');
clims = [ - 1.5 1.5 ];
imagesc( L3H_4, clims )
axis equal
axis( [ 0.5 64.5 0.5 16.5 ] )
title( 'Interpolation Matrix, L_{kl}^{3}(x_{H}^{3})', 'fontsize', 7 )
xlabel( 'N = 64', 'fontsize', 8 )
ylabel( 'N^{3} = 16', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

print( fig_2, '-dpng', 'InterpolationMatricesDims.png' );

% --- Mass and Stiffness Matrices ---

fig_3 = figure( 3 );

subplot( 2, 2, 1 )

Mij_4 = ReadMatrix('Mij_4.dat');
clims = [ -0.035 0.035 ];
imagesc( Mij_4, clims )
axis equal
axis( [ 0.5 64.5 0.5 64.5 ] )
title( 'Mass Matrix, M_{jk}', 'fontsize', 7 )
xlabel( 'N = 64', 'fontsize', 8 )
ylabel( 'N = 64', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap gray

subplot( 2, 2, 2 )

S1ij_4 = ReadMatrix('S1ij_4.dat');
clims = [ -0.2 +0.2 ];
imagesc( S1ij_4, clims )
axis equal
axis( [ 0.5 64.5 0.5 64.5 ] )
title( 'Stiffness Matrix, S_{jk}^{1}', 'fontsize', 7 )
xlabel( 'N = 64', 'fontsize', 8 )
ylabel( 'N = 64', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap gray

subplot( 2, 2, 3 )

S2ij_4 = ReadMatrix('S2ij_4.dat');
clims = [ -0.2 +0.2 ];
imagesc( S2ij_4, clims )
axis equal
axis( [ 0.5 64.5 0.5 64.5 ] )
title( 'Stiffness Matrix, S_{jk}^{2}', 'fontsize', 7 )
xlabel( 'N = 64', 'fontsize', 8 )
ylabel( 'N = 64', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap gray

subplot( 2, 2, 4 )

S3ij_4 = ReadMatrix('S3ij_4.dat');
clims = [ -0.2 +0.2 ];
imagesc( S3ij_4, clims )
axis equal
axis( [ 0.5 64.5 0.5 64.5 ] )
title( 'Stiffness Matrix, S_{jk}^{3}', 'fontsize', 7 )
xlabel( 'N = 64', 'fontsize', 8 )
ylabel( 'N = 64', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap gray

print( fig_3, '-dpng', 'MassAndStiffnessMatrices.png' );

% --- Integration Matrices ---

fig_4 = figure( 4 );

subplot( 2, 3, 1 )

M1ijL_4 = ReadMatrix('M1ijL_4.dat');
clims = [ -0.165 0.165 ];
imagesc( M1ijL_4, clims )
axis equal
axis( [ 0.5 16.5 0.5 64.5 ] )
title( 'Surface Mass Matrix, M_{jk}^{1}(x_{L}^{1})', 'fontsize', 7 )
xlabel( 'N = 16', 'fontsize', 8 )
ylabel( 'N = 64', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap gray

subplot( 2, 3, 2 )

M2ijL_4 = ReadMatrix('M2ijL_4.dat');
clims = [ -0.165 0.165 ];
imagesc( M2ijL_4, clims )
axis equal
axis( [ 0.5 16.5 0.5 64.5 ] )
title( 'Surface Mass Matrix, M_{jk}^{2}(x_{L}^{2})', 'fontsize', 7 )
xlabel( 'N = 16', 'fontsize', 8 )
ylabel( 'N = 64', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap gray

subplot( 2, 3, 3 )

M3ijL_4 = ReadMatrix('M3ijL_4.dat');
clims = [ -0.165 0.165 ];
imagesc( M3ijL_4, clims )
axis equal
axis( [ 0.5 16.5 0.5 64.5 ] )
title( 'Surface Mass Matrix, M_{jk}^{3}(x_{L}^{3})', 'fontsize', 7 )
xlabel( 'N = 16', 'fontsize', 8 )
ylabel( 'N = 64', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap gray

subplot( 2, 3, 4 )

M1ijH_4 = ReadMatrix('M1ijH_4.dat');
clims = [ -0.165 0.165 ];
imagesc( M1ijH_4, clims )
axis equal
axis( [ 0.5 16.5 0.5 64.5 ] )
title( 'Surface Mass Matrix, M_{jk}^{1}(x_{H}^{1})', 'fontsize', 7 )
xlabel( 'N = 16', 'fontsize', 8 )
ylabel( 'N = 64', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap gray

subplot( 2, 3, 5 )

M2ijH_4 = ReadMatrix('M2ijH_4.dat');
clims = [ -0.165 0.165 ];
imagesc( M2ijH_4, clims )
axis equal
axis( [ 0.5 16.5 0.5 64.5 ] )
title( 'Surface Mass Matrix, M_{jk}^{2}(x_{H}^{2})', 'fontsize', 7 )
xlabel( 'N = 16', 'fontsize', 8 )
ylabel( 'N = 64', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap gray

subplot( 2, 3, 6 )

M3ijH_4 = ReadMatrix('M3ijH_4.dat');
clims = [ -0.165 0.165 ];
imagesc( M3ijH_4, clims )
axis equal
axis( [ 0.5 16.5 0.5 64.5 ] )
title( 'Surface Mass Matrix, M_{jk}^{3}(x_{H}^{3})', 'fontsize', 7 )
xlabel( 'N = 16', 'fontsize', 8 )
ylabel( 'N = 64', 'fontsize', 8 )
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
colormap gray

print( fig_4, '-dpng', 'IntegrationMatrices.png' );

% --- Mass Matrices ---

fig_5 = figure( 5 );

Mij_4 = ReadMatrix( 'M0ij_4.dat' );

clims = [ -8 -2 ];
imagesc( log10( abs( Mij_4 ) ), clims )
axis square
title( 'Mass Matrix ( N = 4 ), log_{10}( | M^{(0)} | )', 'fontsize', 15 )
xlabel( 'N = 64', 'fontsize', 20 )
ylabel( 'N = 64', 'fontsize', 20 )
colorbar( 'eastoutside' )
colormap gray
print( fig_5, '-dpng', 'MassMatrix0_G4.png' );

fig_6 = figure( 6 );

Mij_5 = ReadMatrix( 'M0ij_5.dat' );

clims = [ -8 -2 ];
imagesc( log10( abs( Mij_5 ) ), clims )
axis square
title( 'Mass Matrix ( N = 5 ), log_{10}( | M^{(0)} | )', 'fontsize', 15 )
xlabel( 'N = 64', 'fontsize', 20 )
ylabel( 'N = 64', 'fontsize', 20 )
colorbar( 'eastoutside' )
colormap gray
print( fig_6, '-dpng', 'MassMatrix0_G5.png' );