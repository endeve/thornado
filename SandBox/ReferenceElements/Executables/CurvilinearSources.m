x1   = ReadVector(  'x1.dat' );
x2   = ReadVector(  'x2.dat' );
x3   = ReadVector(  'x3.dat' );
J    = ReadVector(   'J.dat' );
g11  = ReadVector( 'g11.dat' );
g22  = ReadVector( 'g22.dat' );
g33  = ReadVector( 'g33.dat' );
sg11 = ReadVector( 'sg11.dat' );
sg22 = ReadVector( 'sg22.dat' );
sg33 = ReadVector( 'sg33.dat' );
dX   = ReadVector(  'dX.dat' );
xL   = ReadVector(  'xL.dat' );
xH   = ReadVector(  'xH.dat' );

L1L = ReadMatrix( 'L1L_4.dat' );
L1H = ReadMatrix( 'L1H_4.dat' );
L2L = ReadMatrix( 'L2L_4.dat' );
L2H = ReadMatrix( 'L2H_4.dat' );
L3L = ReadMatrix( 'L3L_4.dat' );
L3H = ReadMatrix( 'L3H_4.dat' );

dL1 = ReadMatrix( 'dLdX1.dat' );
dL2 = ReadMatrix( 'dLdX2.dat' );
dL3 = ReadMatrix( 'dLdX3.dat' );

M   = ReadMatrix(   'Mij_4.dat' );
M1L = ReadMatrix( 'M1ijL_4.dat' );
M1H = ReadMatrix( 'M1ijH_4.dat' );
M2L = ReadMatrix( 'M2ijL_4.dat' );
M2H = ReadMatrix( 'M2ijH_4.dat' );
M3L = ReadMatrix( 'M3ijL_4.dat' );
M3H = ReadMatrix( 'M3ijH_4.dat' );

S1  = ReadMatrix( 'S1ij_4.dat' );
S2  = ReadMatrix( 'S2ij_4.dat' );
S3  = ReadMatrix( 'S3ij_4.dat' );

[ max(abs(g11-1.0)) ...
  max(abs(g22-(x1).^2)) ...
  max(abs(g33-(x1.*sin(x2)).^2)) ...
  max(abs(J-x1.^2.*sin(x2))) ]'

% --- Compute Volume ---

w  = ReadVector(  'w.dat' );

V_N = prod(dX) * w' * J;
V_A = (1.0/3.0)*(xH(1)^3-xL(1)^3)*(cos(xL(2))-cos(xH(2)))*(xH(3)-xL(3));

[ abs(V_N-V_A) ]'

% --- Compute Areas ---

w1 = ReadVector( 'w1.dat' );
w2 = ReadVector( 'w2.dat' );
w3 = ReadVector( 'w3.dat' );

A1L_N = dX(2)*dX(3)*( w1' * ( L1L * sqrt(g22.*g33) ) );
A1H_N = dX(2)*dX(3)*( w1' * ( L1H * sqrt(g22.*g33) ) );
A1L_A = xL(1)^2*(cos(xL(2))-cos(xH(2)))*(xH(3)-xL(3));
A1H_A = xH(1)^2*(cos(xL(2))-cos(xH(2)))*(xH(3)-xL(3));

[ abs(A1L_N-A1L_A) abs(A1H_N-A1H_A) ]'

A2L_N = dX(1)*dX(3)*( w2' * ( L2L * sqrt(g11.*g33) ) );
A2H_N = dX(1)*dX(3)*( w2' * ( L2H * sqrt(g11.*g33) ) );

A2L_A = (1.0/2.0)*(xH(1)^2-xL(1)^2)*sin(xL(2))*(xH(3)-xL(3));
A2H_A = (1.0/2.0)*(xH(1)^2-xL(1)^2)*sin(xH(2))*(xH(3)-xL(3));

[ abs(A2L_N-A2L_A) abs(A2H_N-A2H_A) ]'

A3L_N = dX(1)*dX(2)*( w3' * ( L3L * sqrt(g11.*g22) ) );
A3H_N = dX(1)*dX(2)*( w3' * ( L3H * sqrt(g11.*g22) ) );

A3L_A = (1.0/2.0)*(xH(1)^2-xL(1)^2)*(xH(2)-xL(2));
A3H_A = (1.0/2.0)*(xH(1)^2-xL(1)^2)*(xH(2)-xL(2));

[ abs(A3L_N-A3L_A) abs(A3H_N-A3H_A) ]'

% --- 

dg22dx1 = M \ ( M1H * ( L1H * sg22 ) - M1L * ( L1L * sg22 ) - S1 * sg22 );
dg33dx1 = M \ ( M1H * ( L1H * sg33 ) - M1L * ( L1L * sg33 ) - S1 * sg33 );

dg22dx2 = M \ ( M2H * ( L2H * sg22 ) - M2L * ( L2L * sg22 ) - S2 * sg22 );
dg33dx2 = M \ ( M2H * ( L2H * sg33 ) - M2L * ( L2L * sg33 ) - S2 * sg33 );

dg22dx3 = M \ ( M3H * ( L3H * sg22 ) - M3L * ( L3L * sg22 ) - S3 * sg22 );
dg33dx3 = M \ ( M3H * ( L3H * sg33 ) - M3L * ( L3L * sg33 ) - S3 * sg33 );

Div1_1 = ( w1' * ( L1H * J ) - w1' * ( L1L * J ) );
Div1_2 = w' * ( (J./sg22).*dg22dx1 + (J./sg33).*dg33dx1 );

Div2_1 = ( w2' * ( L2H * J ) - w2' * ( L2L * J ) );
Div2_2 = w' * ( (J./sg22).*dg22dx2 + (J./sg33).*dg33dx2 );

Div3_1 = ( w3' * ( L3H * J ) - w3' * ( L3L * J ) );
Div3_2 = w' * ( (J./sg22).*dg22dx3 + (J./sg33).*dg33dx3 );