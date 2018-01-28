% close all
clear all

nPts = 1000000;
A = zeros(nPts,1);
B = zeros(nPts,1);
for iPt = 1 : nPts
  
  A(iPt) = 1.d2*sign(2*(0.5-rand))*10^(6*(rand-1.0));
  B(iPt) = 1.d2*sign(2*(0.5-rand))*10^(6*(rand-1.0));
  
end

nMu = 256;
MuL = - 0.999999;
MuR = + 0.999999;
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

tic
for iPt = 1 : nPts
  
  J_FD(iPt)     = + 0.5 * Trapez( Mu, Mu.^0 .* f_FD(:,iPt) );
  H_FD(iPt)     = + 0.5 * Trapez( Mu, Mu.^1 .* f_FD(:,iPt) );
 
end
toc

clf
plot( H_FD, J_FD, '.k', 'MarkerSize', 6 );hold on
yy = linspace(0.0,1.0,1024)';
xxP_FD = + yy.*(1.0-yy);
xxM_FD = - yy.*(1.0-yy);
xxP_MB = + yy;
xxM_MB = - yy;
plot( xxP_MB, yy, '-r', 'linewidth', 3 )
plot( xxM_MB, yy, '-r', 'linewidth', 3 )
plot( xxP_FD, yy, '-b', 'linewidth', 3 )
plot( xxM_FD, yy, '-b', 'linewidth', 3 )
hold off
axis([-0.5 0.5 0.0 1.0]);
xlabel('Angular Moment H');
ylabel('Angular Moment J');
title('Realizble Set R');
