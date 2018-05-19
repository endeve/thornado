close all
clear all

fig = figure;
yy = linspace(0.0,1.0,1024)'; hold on
xxP_FD = + yy.*(1.0-yy);
xxM_FD = - yy.*(1.0-yy);
xxP_MB = + yy;
xxM_MB = - yy;
fill([xxP_FD flip(xxP_FD)],[yy zeros(size(yy))],[0.7 0.8 1.0],'LineStyle','none')
fill([xxM_FD flip(xxM_FD)],[yy zeros(size(yy))],[0.7 0.8 1.0],'LineStyle','none')
plot( xxP_MB, yy, '-r', 'linewidth', 2 )
plot( xxM_MB, yy, '-r', 'linewidth', 2 )
plot( xxP_FD, yy, '-k', 'linewidth', 2 )
plot( xxM_FD, yy, '-k', 'linewidth', 2 )
hold off
axis([-0.5 0.5 0.0 1.0]);
xlabel('Angular Moment H');
ylabel('Angular Moment J');
title('Realizble Set R');
box on
saveas(fig,'RealizableSetFermionic.png');
return