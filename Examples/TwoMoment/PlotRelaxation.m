close all
clear all

iFile = 11;

[ t, E, X1, X2, X3, N_1 ] = ReadRadiationFields( 'Relaxation', iFile, 1 );
[ ~, ~,  ~,  ~,  ~, N_2 ] = ReadRadiationFields( 'Relaxation', iFile, 2 );
[ ~, ~,  ~,  ~,  ~, T, ~, ~, ~, Me, Mp, Mn ] = ReadFluidFields_Auxiliary( 'Relaxation', iFile );

kT    = 8.617333262d-11 * T(1,1,1);
Mnu_1 = Me(1,1,1) + Mp(1,1,1) - Mn(1,1,1);
Mnu_2 = - Mnu_1;

N0_1 = 1.0 ./ ( exp( ( E - Mnu_1 ) ./ kT ) + 1.0 );
N0_2 = 1.0 ./ ( exp( ( E - Mnu_2 ) ./ kT ) + 1.0 );

fig_1 = figure( 1 );

subplot( 2, 1, 1 );

loglog( E, N0_1        , '--b', 'linewidth', 2 ); hold on
loglog( E, N_1(:,1,1,1),  'ob', 'linewidth', 2 )
loglog( E, N0_2        , '--r', 'linewidth', 2 )
loglog( E, N_2(:,1,1,1),  'or', 'linewidth', 2 )

subplot( 2, 1, 2);

loglog( E, abs(N0_1-N_1(:,1,1,1))./N0_1, '-b', 'linewidth', 2 ); hold on
loglog( E, abs(N0_2-N_2(:,1,1,1))./N0_2, '-r', 'linewidth', 2 )