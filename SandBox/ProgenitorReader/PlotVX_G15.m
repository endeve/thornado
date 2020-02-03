close all
clear all

R_003 = ReadVector( './Executables/R_003ms.dat' );
R_010 = ReadVector( './Executables/R_010ms.dat' );
R_050 = ReadVector( './Executables/R_050ms.dat' );
R_100 = ReadVector( './Executables/R_100ms.dat' );
R_150 = ReadVector( './Executables/R_150ms.dat' );
R_250 = ReadVector( './Executables/R_250ms.dat' );

D_003 = ReadVector( './Executables/D_003ms.dat' );
D_010 = ReadVector( './Executables/D_010ms.dat' );
D_050 = ReadVector( './Executables/D_050ms.dat' );
D_100 = ReadVector( './Executables/D_100ms.dat' );
D_150 = ReadVector( './Executables/D_150ms.dat' );
D_250 = ReadVector( './Executables/D_250ms.dat' );

fig_1 = figure( 1 );

semilogy( R_003/1.d5, D_003, '-k', 'linewidth', 3 ); hold on
semilogy( R_010/1.d5, D_010, '-b', 'linewidth', 3 )
semilogy( R_050/1.d5, D_050, '-r', 'linewidth', 3 )
semilogy( R_100/1.d5, D_100, '-m', 'linewidth', 3 )
semilogy( R_150/1.d5, D_150, '-c', 'linewidth', 3 )
semilogy( R_250/1.d5, D_250, '-g', 'linewidth', 3 )
set( gca, 'fontsize',   20 )
set( gca, 'linewidth',  02 )
xlabel( 'Radius [ km ]', 'fontsize', 20 )
ylabel( 'Mass Density [ g cm^{-3} ]', 'fontsize', 20 )
axis( [ 0 200 1.d7 1.d15 ] )
legend_h...
  = legend...
      ( '003 ms', '010 ms', '050 ms', '100 ms', '150 ms', '250 ms',...
        'location', 'northeast' );
set( legend_h, 'fontsize', 20 );

print( fig_1, '-dpng', 'MassDensity_VX_G15.png' )

T_003 = ReadVector( './Executables/T_003ms.dat' );
T_010 = ReadVector( './Executables/T_010ms.dat' );
T_050 = ReadVector( './Executables/T_050ms.dat' );
T_100 = ReadVector( './Executables/T_100ms.dat' );
T_150 = ReadVector( './Executables/T_150ms.dat' );
T_250 = ReadVector( './Executables/T_250ms.dat' );

fig_2 = figure( 2 );

semilogy( R_003/1.d5, T_003, '-k', 'linewidth', 3 ); hold on
semilogy( R_010/1.d5, T_010, '-b', 'linewidth', 3 )
semilogy( R_050/1.d5, T_050, '-r', 'linewidth', 3 )
semilogy( R_100/1.d5, T_100, '-m', 'linewidth', 3 )
semilogy( R_150/1.d5, T_150, '-c', 'linewidth', 3 )
semilogy( R_250/1.d5, T_250, '-g', 'linewidth', 3 )
set( gca, 'fontsize',   20 )
set( gca, 'linewidth',  02 )
xlabel( 'Radius [ km ]', 'fontsize', 20 )
ylabel( 'Temperature [ MeV ]', 'fontsize', 20 )
axis( [ 0 200 1.d-1 1.d+2 ] )
legend_h...
  = legend...
      ( '003 ms', '010 ms', '050 ms', '100 ms', '150 ms', '250 ms',...
        'location', 'southwest' );
set( legend_h, 'fontsize', 20 );

print( fig_2, '-dpng', 'Temperature_VX_G15.png' )

Y_003 = ReadVector( './Executables/Y_003ms.dat' );
Y_010 = ReadVector( './Executables/Y_010ms.dat' );
Y_050 = ReadVector( './Executables/Y_050ms.dat' );
Y_100 = ReadVector( './Executables/Y_100ms.dat' );
Y_150 = ReadVector( './Executables/Y_150ms.dat' );
Y_250 = ReadVector( './Executables/Y_250ms.dat' );

fig_3 = figure( 3 );

plot( R_003/1.d5, Y_003, '-k', 'linewidth', 3 ); hold on
plot( R_010/1.d5, Y_010, '-b', 'linewidth', 3 )
plot( R_050/1.d5, Y_050, '-r', 'linewidth', 3 )
plot( R_100/1.d5, Y_100, '-m', 'linewidth', 3 )
plot( R_150/1.d5, Y_150, '-c', 'linewidth', 3 )
plot( R_250/1.d5, Y_250, '-g', 'linewidth', 3 )
set( gca, 'fontsize',   20 )
set( gca, 'linewidth',  02 )
xlabel( 'Radius [ km ]', 'fontsize', 20 )
ylabel( 'Electron Fraction', 'fontsize', 20 )
axis( [ 0 200 0.05 0.6 ] )
legend_h...
  = legend...
      ( '003 ms', '010 ms', '050 ms', '100 ms', '150 ms', '250 ms',...
        'location', 'southeast' );
set( legend_h, 'fontsize', 20 );

print( fig_3, '-dpng', 'ElectronFraction_VX_G15.png' )