close all
clear all

PolynomialCB = @(x) x^2 * ( 3.0 - abs(x) + 3.0 * x^2 ) / 5.0;
EddingtonFCB = @(e,f) (2.0/3.0)*(1.0-e)*(1.0-2.0*e)*PolynomialCB(f/(1.0-e))+1.0/3.0;
PolynomialBL = @(x) ( 9.0 * x^2 - 5.0 + sqrt( 33.0 * x^4 - 42.0 * x^2 + 25.0 ) ) /8.0;
EddingtonFBL = @(e,f) (2.0/3.0)*(1.0-e)*(1.0-2.0*e)*PolynomialBL(f/(1.0-e))+1.0/3.0;
PolynomialKS = @(x) x^2;
EddingtonFKS = @(e,f) (2.0/3.0)*(1.0-e)*(1.0-2.0*e)*PolynomialKS(f/(1.0-e))+1.0/3.0; 

nPts = 2000000;
J_FD = zeros(nPts,1);
H_FD = zeros(nPts,1);
iPt = 1;
while( iPt <= nPts )
    J_FD(iPt) = rand(1);
    temp = (rand(1)-0.5)*0.5;
    if( abs(temp) <= (1-J_FD(iPt))*J_FD(iPt) )
        H_FD(iPt) = temp;
        iPt = iPt + 1;
    end
end

K_CB = zeros(nPts,1);
K_BL = zeros(nPts,1);
K_MI = zeros(nPts,1);
K_KS = zeros(nPts,1);

tic
for iPt = 1 : nPts
  
  K_CB(iPt) = EddingtonFCB( J_FD(iPt), H_FD(iPt)/J_FD(iPt) ) * J_FD(iPt);
  K_BL(iPt) = EddingtonFBL( J_FD(iPt), H_FD(iPt)/J_FD(iPt) ) * J_FD(iPt);
  K_MI(iPt) = EddingtonFCB( 0.0,       H_FD(iPt)/J_FD(iPt) ) * J_FD(iPt);
  K_KS(iPt) = EddingtonFBL( J_FD(iPt), H_FD(iPt)/J_FD(iPt) ) * J_FD(iPt);
  
end
toc

ClosureName = {'Kershaw', 'CB', 'BL','Minerbo'};
ClosureSaveName = {'BLKS','CBME','BLME','MI'};

for ii = 1:max(size(ClosureName))
    
Closure = ClosureName{ii}; % 
FileName = ClosureSaveName{ii};
switch Closure
    case 'Minerbo'
        K_plot = K_MI;
    case 'CB'
        K_plot = K_CB;
    case 'BL'
        K_plot = K_BL;
    case 'Kershaw'
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
axis([-0.3 0.300001 0.0-1.d-2 1.0+1.d-2]);
set(gca,'FontSize',15);
xlabel('$\mathcal{H}_{ab}$','Interpreter','LaTeX','Fontsize',19);
ylabel('$\mathcal{J}_{ab}$','Interpreter','LaTeX','Fontsize',19);
title(['\bf ',Closure],'Interpreter','LaTeX','Fontsize',19);
x0=10;
y0=5;
width=400;
height=320;
set(gcf,'units','points','position',[x0,y0,width,height])
pause(1)
saveas(fig,['MabWith', FileName,'.png']);
end