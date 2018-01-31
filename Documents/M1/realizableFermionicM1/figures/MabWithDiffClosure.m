% close all
clear all

nPts = 1000000;


PolynomialCB = @(x) x^2 * ( 3.0 - abs(x) + 3.0 * x^2 ) / 5.0;
EddingtonFCB = @(e,f) (2.0/3.0)*(1.0-e)*(1.0-2.0*e)*PolynomialCB(f/(1.0-e))+1.0/3.0;
PolynomialBL = @(x) ( 9.0 * x^2 - 5.0 + sqrt( 33.0 * x^4 - 42.0 * x^2 + 25.0 ) ) /8.0;
EddingtonFBL = @(e,f) (2.0/3.0)*(1.0-e)*(1.0-2.0*e)*PolynomialBL(f/(1.0-e))+1.0/3.0;
PolynomialKS = @(x) x^2;
EddingtonFKS = @(e,f) (2.0/3.0)*(1.0-e)*(1.0-2.0*e)*PolynomialKS(f/(1.0-e))+1.0/3.0; 

A = zeros(nPts,1);
B = zeros(nPts,1);
for iPt = 1 : nPts
  
  A(iPt) = 1.d2*sign(2*(0.5-rand))*10^(6*(rand-1.0));
  B(iPt) = 1.d2*sign(2*(0.5-rand))*10^(6*(rand-1.0));
  
end

nMu = 256;
MuL = - 1.0;
MuR = + 1.0;
Mu  = linspace( MuL, MuR, nMu )';
dMu = ( MuR - MuL ) / nMu;

% Compute Angular Distributions for Various Parameters:

f_FD = zeros(nMu,nPts);

for iPt = 1 : nPts
  for iMu = 1 : nMu
    
    f_FD(iMu,iPt) = 1.0/(exp(A(iPt)+B(iPt)*Mu(iMu))+1.0); % Fermi-Dirac
    
  end
end

% Compute Angular Moments:

J_FD     = zeros(nPts,1);
H_FD     = zeros(nPts,1);
K_FD     = zeros(nPts,1);
K_CB = zeros(nPts,1);
K_BL = zeros(nPts,1);
K_MI = zeros(nPts,1);
K_KS = zeros(nPts,1);


tic
for iPt = 1 : nPts
  
  J_FD(iPt)     = + 0.5 * Trapez( Mu, Mu.^0 .* f_FD(:,iPt) );
  H_FD(iPt)     = + 0.5 * Trapez( Mu, Mu.^1 .* f_FD(:,iPt) );
  K_FD(iPt)     = + 0.5 * Trapez( Mu, Mu.^2 .* f_FD(:,iPt) );
  K_CB(iPt) = EddingtonFCB( J_FD(iPt), H_FD(iPt)/J_FD(iPt) ) * J_FD(iPt);
  K_BL(iPt) = EddingtonFBL( J_FD(iPt), H_FD(iPt)/J_FD(iPt) ) * J_FD(iPt);
  K_MI(iPt) = EddingtonFCB( 0.0,       H_FD(iPt)/J_FD(iPt) ) * J_FD(iPt);
  K_KS(iPt) = EddingtonFBL( J_FD(iPt), H_FD(iPt)/J_FD(iPt) ) * J_FD(iPt);
  
end
toc

f_FD = 0;
f_MB = 0;

ClosureName = {'BLKS', 'CBME', 'BLME','MI'};

for ii = 1:max(size(ClosureName))
    
Closure = ClosureName{ii}; % 

switch Closure
    case 'MI'
        K_plot = K_MI;
    case 'CBME'
        K_plot = K_CB;
    case 'BLME'
        K_plot = K_BL;
    case 'BLKS'
        K_plot = K_KS;
end

M1 = 0.5 .* ( J_FD(1:end-1) +   H_FD(1:end-1) ...
            + J_FD(2:end)   -   H_FD(2:end) );
M2 = 0.5 .* ( H_FD(1:end-1) + K_plot(1:end-1)...
            + H_FD(2:end)   - K_plot(2:end) );

fig = figure;
plot( M2, M1, '.','Color',[0.7 0.8 1.0], 'MarkerSize', 5.5 );hold on
yy = linspace(0.0,1.0,1024)';
xxP_FD = + yy.*(1.0-yy);
xxM_FD = - yy.*(1.0-yy);
xxP_MB = + yy;
xxM_MB = - yy;
plot( xxP_FD, yy, '-k', 'linewidth', 2 )
plot( xxM_FD, yy, '-k', 'linewidth', 2 )
hold off
axis([-0.5 0.5 0.0 1.0]);
xlabel('H_{ab}');
ylabel('J_{ab}');
title(['M_{ab} with ', Closure,' closure']);
pause(1)
saveas(fig,['MabWith', Closure,'.png']);
end