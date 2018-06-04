close all
clear all

fig = figure;
yy = linspace(0.0,1.0,1024)'; hold on
xxP_FD = + yy.*(1.0-yy);
xxM_FD = - yy.*(1.0-yy);

xxP_MB = + yy;
xxM_MB = - yy;
xxP_MB(end) = 0;
xxM_MB(end) = 0;
fill([xxP_MB xxM_MB], [yy yy], [1.0 0.9 0.9],'LineStyle','none')
alpha(0.8)

xxP_MB = + yy;
xxM_MB = - yy;
fill([xxP_FD flip(xxM_FD)],[yy yy],[0.7 0.8 1.0],'LineStyle','none')
plot( xxP_MB, yy, '-r', 'linewidth', 2 )
plot( xxM_MB, yy, '-r', 'linewidth', 2 )
plot( xxP_FD, yy, '-k', 'linewidth', 2 )
plot( xxM_FD, yy, '-k', 'linewidth', 2 )
hold off
axis([-0.5 0.5 0.0 1.0]);
set(gca,'FontSize',10)
set(gca,'xtick',-0.5:0.1:0.5);
xlabel('{\bf Angular Moment} $\mathcal{H}$','Interpreter','LaTeX','Fontsize',11);
ylabel('{\bf Angular Moment} $\mathcal{J}$','Interpreter','LaTeX','Fontsize',11);
title('{\bf Realizable set} $\mathcal{R}$','Interpreter','LaTeX','Fontsize',13)
box on
set(gca, 'Layer', 'top')
saveas(fig,'RealizableSetFermionic.png');
return
