close all
clear all

% Spatial Grid:
minR  = 0;
maxR  = 2;
nR    = 256;
R     = logspace( minR, maxR, nR )';

% Profile Parameters:
R_D = 20.0;
H_D = 10.0;

R_T = 25.0;
H_T = 20.0;

R_Y = 45.0;
H_Y = 10.0;

% Mass Density:
minD = 1.0d08; % g / cm^3
maxD = 4.0d14; % g / cm^3
D    = 0.5.*(maxD.*(1.0-tanh((R-R_D)./H_D))+minD.*(1.0-tanh((R_D-R)./H_D)))';

% Temperature:
minT = 5.0d09; % K
maxT = 2.6d11; % K
T    = 0.5.*(maxT.*(1.0-tanh((R-R_T)./H_T))+minT.*(1.0-tanh((R_T-R)./H_T)))';

% Electron Fraction:
minY = 0.30;
maxY = 0.46;
Y    = 0.5.*(maxY.*(1.0+tanh((R-R_Y)./H_Y))+minY.*(1.0+tanh((R_Y-R)./H_Y)))';

[ D_T, T_T, Y_T, nD_T, nT_T, nY_T, P_T, S_T, E_T, Me_T, Mp_T, Mn_T,...
  Xp_T, Xn_T, Xa_T, Xh_T, Zh_T, Ah_T, Eh_T, Eth_T, Gm_T, OS_EOS ]...
  = readEosTable( 'wl-EOS-SFHo-15-25-50.h5' );

P  = zeros( nR, 1 );
S  = zeros( nR, 1 );
E  = zeros( nR, 1 );
Me = zeros( nR, 1 );
Mp = zeros( nR, 1 );
Mn = zeros( nR, 1 );
Xp = zeros( nR, 1 );
Xn = zeros( nR, 1 );
Xa = zeros( nR, 1 );
Xh = zeros( nR, 1 );
for i = 1 : nR

  P(i)  = interpolateEos( D(i), T(i), Y(i), D_T, T_T, Y_T, P_T,  OS_EOS(01) );
  S(i)  = interpolateEos( D(i), T(i), Y(i), D_T, T_T, Y_T, S_T,  OS_EOS(02) );
  E(i)  = interpolateEos( D(i), T(i), Y(i), D_T, T_T, Y_T, E_T,  OS_EOS(03) );
  Me(i) = interpolateEos( D(i), T(i), Y(i), D_T, T_T, Y_T, Me_T, OS_EOS(04) );
  Mp(i) = interpolateEos( D(i), T(i), Y(i), D_T, T_T, Y_T, Mp_T, OS_EOS(05) );
  Mn(i) = interpolateEos( D(i), T(i), Y(i), D_T, T_T, Y_T, Mn_T, OS_EOS(06) );
  Xp(i) = interpolateEos( D(i), T(i), Y(i), D_T, T_T, Y_T, Xp_T, OS_EOS(07) );
  Xn(i) = interpolateEos( D(i), T(i), Y(i), D_T, T_T, Y_T, Xn_T, OS_EOS(08) );
  Xa(i) = interpolateEos( D(i), T(i), Y(i), D_T, T_T, Y_T, Xa_T, OS_EOS(09) );
  Xh(i) = interpolateEos( D(i), T(i), Y(i), D_T, T_T, Y_T, Xh_T, OS_EOS(10) );
  
end

fig_1 = figure(1);
semilogy( R, D, '-k', 'linewidth', 2 )
axis( [ 0.0d0 1.0d2 1.0d08 1.0d15 ] );
xlabel( 'Radius [ km ]', 'fontsize', 15 )
ylabel( 'Mass Density [ g cm^{-3} ]', 'fontsize', 15 )
print( fig_1, '-dpng', [ './MockProfile_Density.png' ] )

fig_2 = figure(2);
semilogy( R, T, '-k', 'linewidth', 2 )
axis( [ 0.0d0 1.0d2 3.0d09 3.0d11 ] );
xlabel( 'Radius [ km ]', 'fontsize', 15 )
ylabel( 'Temperature [ K ]', 'fontsize', 15 )
print( fig_2, '-dpng', [ './MockProfile_Temperature.png' ] )

fig_3 = figure(3);
plot( R, Y, '-k', 'linewidth', 2 )
axis( [ 0.0d0 1.0d2 0.25 0.5 ] );
xlabel( 'Radius [ km ]', 'fontsize', 15 )
ylabel( 'Electron Fraction', 'fontsize', 15 )
print( fig_3, '-dpng', [ './MockProfile_ElectronFraction.png' ] )

fig_4 = figure(4);
plot...
  ( R, Xp,  '-k',...
    R, Xn,  '-b',...
    R, Xa,  ':r',...
    R, Xh, '--g',...
    R, Xp+Xn+Xa+Xh, '--k',...
    'linewidth',2 )
axis( [ 0.0d0 1.0d2 -0.05 1.05 ] );
xlabel( 'Radius [ km ]', 'fontsize', 15 )
ylabel( 'Mass Fractions', 'fontsize', 15 )
legend_h...
  = legend( 'X_p', 'X_n', 'X_{\alpha}', 'X_h',...
            'location', 'east' );
set( legend_h, 'fontsize', 20 );
print( fig_4, '-dpng', [ './MockProfile_MassFractions.png' ] )

%
% Neutrino Opacities:

% Energy Grid:
minEnu = 0;
maxEnu = log10(2.9d2);
nEnu   = 40;
Enu    = logspace( minEnu, maxEnu, nEnu )';

[ E_T, D_T, T_T, Y_T, Chi_T, OS_Chi ]...
  = readAbsorptionOpacityTable( 'wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5' );
[ E_T, D_T, T_T, Y_T, Sig_T0, Sig_T1, OS_Sig ]...
  = readIsoScatteringOpacityTable( 'wl-Op-SFHo-15-25-50-E40-B85-Iso.h5' );

Chi = zeros( nEnu, nR );
Sig = zeros( nEnu, nR );
f_0 = zeros( nEnu, nR );
kB  = 1.38065d-23;
MeV = 1.60218d-13;
for j = 1 : nR
  Mnu = Mp(j)+Me(j)-Mn(j);
  kT  = kB * T(j) / MeV;
  for i = 1 : nEnu
  
    Chi(i,j)...
      = interpolate1D3D...
          ( Enu(i), D(j), T(j), Y(j), E_T, D_T, T_T, Y_T, [1 1 1 0]', Chi_T, OS_Chi );
    
    Sig(i,j)...
      = interpolate1D3D...
          ( Enu(i), D(j), T(j), Y(j), E_T, D_T, T_T, Y_T, [1 1 1 0]', Sig_T0, OS_Sig(1) );
    
    f_0(i,j)...
      = 1.0 / ( exp( ( Enu(i) - Mnu ) / kT ) + 1.0 );
  
  end
end

Chi_A = zeros( nR, 1 );
Sig_A = zeros( nR, 1 );
for i = 1 : nR

  Chi_A(i)...
    = Trapez( Enu(:), Chi(:,i) .* f_0(:,i) .* Enu(:).^2 )...
        / Trapez( Enu(:), f_0(:,i) .* Enu(:).^2 );
  
  Sig_A(i)...
    = Trapez( Enu(:), Sig(:,i) .* f_0(:,i) .* Enu(:).^2 )...
        / Trapez( Enu(:), f_0(:,i) .* Enu(:).^2 );
  
end

Tau = zeros( nR, 1 );
for i = 1 : nR
  Tau(i) = Trapez( R(i:nR), Chi_A(i:nR)+Sig_A(i:nR) )*1.d5;
end

fig_5 = figure(5);
semilogy...
  ( R, 1.0d-5./Chi_A, '-k',...
    R, 1.0d-5./Sig_A, '-b',...
    'linewidth', 2 )
axis( [ 0.0d0 1.0d2 1.d-5 1.d10 ] );
xlabel( 'Radius [ km ]', 'fontsize', 15 )
ylabel( 'Mean Free Path [km]', 'fontsize', 15 )
legend_h...
  = legend( '\lambda_{\chi}', '\lambda_{\sigma}',...
            'location', 'northwest' );
set( legend_h, 'fontsize', 20 );
print( fig_5, '-dpng', [ './MockProfile_MeanFreePath.png' ] )

fig_6 = figure(6);
ind = find(Tau<1.,1,'first');
semilogy...
  ( R, Tau, '-k',...
    [ R(ind) R(ind) ],[1.d-6 1.d+6], '--k',...
    'linewidth', 2 )
axis( [ 0.0d0 1.0d2 1.d-6 1.d+6 ] );
xlabel( 'Radius [ km ]', 'fontsize', 15 )
ylabel( 'Optical Depth', 'fontsize', 15 )