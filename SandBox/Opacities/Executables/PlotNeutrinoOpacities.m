close all
clear all

E = ReadVector( 'E.dat' );

Chi_NuE      = ReadVector(      'Chi_NuE.dat' );
Sig_NuE      = ReadVector(    'Sigma_NuE.dat' );
Chi_NES_NuE  = ReadVector(  'Chi_NES_NuE.dat' );
Chi_Pair_NuE = ReadVector( 'Chi_Pair_NuE.dat' );

Chi_NuE_Bar      = ReadVector(      'Chi_NuE_Bar.dat' );
Sig_NuE_Bar      = ReadVector(    'Sigma_NuE_Bar.dat' );
Chi_NES_NuE_Bar  = ReadVector(  'Chi_NES_NuE_Bar.dat' );
Chi_Pair_NuE_Bar = ReadVector( 'Chi_Pair_NuE_Bar.dat' );

loglog( E, Chi_NuE,              '-k', 'linewidth', 2 ); hold on
loglog( E, Sig_NuE,              '-b', 'linewidth', 2 )
loglog( E, Chi_NES_NuE,          '-c', 'linewidth', 2 )
loglog( E, 1.d3 .* Chi_Pair_NuE, '-r', 'linewidth', 2 )

loglog( E, Chi_NuE_Bar,              '--k', 'linewidth', 3 )
loglog( E, Sig_NuE_Bar,              '--b', 'linewidth', 3 )
loglog( E, Chi_NES_NuE_Bar,          '--c', 'linewidth', 3 )
loglog( E, 1.d3 .* Chi_Pair_NuE_Bar, '--r', 'linewidth', 3 )

set( gca, 'fontsize', 13 )
set( gca, 'fontweight', 'bold' )
set( gca, 'linewidth', 02 )
axis( [ 0.95 300 1e-6 1e-0 ] )
title( 'Neutrino Opacities' )
xlabel( '[ MeV ]',     'fontsize', 18 )
ylabel( '[ cm^{-1} ]', 'fontsize', 18 )