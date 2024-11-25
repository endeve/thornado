clear all
close all

theta_Kershaw = @(x) x.*x;
theta_MECB = @(x) x.*x.*(3-x+3.*x.*x)./5;
theta_MEBL = @(x) (9.*x.*x-5+sqrt(33.*x.^4 - 42.*x.*x+25))/8;

nPts = 1024;
J = 0.01;
h = linspace(0.0,1.0 - J,nPts);

Kershaw = zeros(nPts,1);
MECB = zeros(nPts,1);
MEBL = zeros(nPts,1);
Minerbo = zeros(nPts,1);
bryup = zeros(nPts,1);
brydown = zeros(nPts,1);

Kershaw = Closure(J,h,theta_Kershaw);
MECB = Closure(J,h,theta_MECB);
MEBL = Closure(J,h,theta_MEBL);
Minerbo = 1/3 + 2.*h.*h.*...
                (3-h+3.*h.*h)./15;
bryup = min(ones(1,nPts),1/(3*J)-J.*h.*h./(1-J) );
brydown = max( (1-2/(3*J)).*ones(1,nPts),h.*h ); 
fig = plot(h,Kershaw,'-','linewidth', 2); hold on
plot(h,MECB,'-','linewidth', 2); 
plot(h,MEBL,'--','linewidth', 2); 
plot(h,Minerbo,'-.','linewidth', 2);
plot(h,bryup,'-k','linewidth', 2);
plot(h,brydown,'-k','linewidth', 2);
lgnd = legend({'Kershaw','CB','BL','Minerbo'},...
    'Location','best');
set(lgnd,'FontSize',14);
% xlim([0,1-J]);
set(gca,'xtick',0:0.2:1.0);
set(gca,'FontSize',14)
xlabel('{\bf Flux factor} $h$','Interpreter','LaTeX','FontSize',19);
ylabel('{\bf Eddington factor} $\chi$','Interpreter','LaTeX','FontSize',19);
title(['$\mathcal{J}$ = ',num2str(J)],'Interpreter','LaTeX','FontSize',20);

x0=10;
y0=5;
width=400;
height=320;
set(gcf,'units','points','position',[x0,y0,width,height])
saveas(fig,['Closures0_0', num2str(J*100),'.png']);

function [ value ] = Closure( J, h, theta_fun )

temp = 2.*(1-J).*(1-2.*J)/3 .* theta_fun(h/(1-J));
value = 1/3 + temp;

end