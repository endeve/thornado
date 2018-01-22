close all
clear all

dLdX1_q = ReadMatrix('./Executables/dLdX1_q.dat');
dLdX2_q = ReadMatrix('./Executables/dLdX2_q.dat');
dLdX3_q = ReadMatrix('./Executables/dLdX3_q.dat');
L_X1_Dn = ReadMatrix('./Executables/L_X1_Dn.dat');
L_X2_Dn = ReadMatrix('./Executables/L_X2_Dn.dat');
L_X3_Dn = ReadMatrix('./Executables/L_X3_Dn.dat');

fig_1 = figure(1);

subplot( 2, 3, 1 )

spy( dLdX1_q, 20 );
axis square
title( 'dL/dX^{1}', 'fontsize', 10 )

subplot( 2, 3, 2 )

spy( dLdX2_q, 20 );
axis square
title( 'dL/dX^{2}', 'fontsize', 10 )

subplot( 2, 3, 3 )

spy( dLdX3_q, 20 );
axis square
title( 'dL/dX^{3}', 'fontsize', 10 )

subplot( 2, 3, 4 )

spy( L_X1_Dn, 20 );
axis square
title( 'L(X^{-},Y,Z)', 'fontsize', 10 )

subplot( 2, 3, 5 )

spy( L_X2_Dn, 20 );
axis square
title( 'L(X,Y^{-},Z)', 'fontsize', 10 )

subplot( 2, 3, 6 )

spy( L_X3_Dn, 20 );
axis square
title( 'L(X,Y,Z^{-})', 'fontsize', 10 )

print...
  ( fig_1, '-dpng', 'ElementMatrices.png' )