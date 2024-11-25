clear all
close all

kappa = @(x) x.*x.*(3-x+3.*x.*x)./5;
Janka = @(a, b, n, h) (1+a.*power(h,b)+(2-a).*power(h,n))./3;

nPts = 20;
Kershaw = zeros(nPts,1);
Wilson = zeros(nPts,1);
Levermore = zeros(nPts,1);
Minerbo = zeros(nPts,1); % ME
CB = zeros(nPts,1);    % MEFD
Janka1 = zeros(nPts,1);
Janka2 = zeros(nPts,1);
data = zeros(nPts,7);

bryup = zeros(nPts,1);
brydown = zeros(nPts,1);

% legendArr = {'Kershaw','Wilson','Levermore','Minerbo',...
%     'CB','Janka 1','Janka 2'};
% marker = {'-','-','-','-','--','-','-'};

J_arr = [0.1, 0.9]; %0.01, 0.4, 0.6, 0.99

load('LegendColorMarker.mat');
color(:,4) = [0.2,0.5,1];

for kk = 1:length(J_arr)
%     close all
    J = J_arr(kk);
    h = linspace(0.0,1.0 - J,nPts);
    
    Kershaw = (1 + 2.*h.*h)./3;
    Wilson = 1/3 - h./3 + h.*h;
    Levermore = (3 + 4.*h.*h)./(5 + 2.*sqrt(4-3.*h.*h));
    Minerbo = 1/3 + 2.*h.*h.*...
        (3-h+3.*h.*h)./15;  
%     MEFD_bry = (1-2.*h+4.*h.*h)./3;
    CB = Closure(J,h,kappa);
    Janka1 = Janka(0.5, 1.3064, 4.1342, h);
    Janka2 = Janka( 1,  1.345,  5.1717, h);
    
    bryup = min(ones(1,nPts),1/(3*J)-J.*h.*h./(1-J) );
    brydown = max( (1-2/(3*J)).*ones(1,nPts),h.*h );
    data = [Kershaw(:),Wilson(:),Levermore(:),Minerbo(:),CB(:),...
        Janka1(:),Janka2(:)];
    fig = figure;
    for ii = 1:7
        plot(h,data(:,ii),marker{ii},'Color',color(:,ii),...
            'linewidth', 2); hold on
    end
    plot(h,bryup,'--k','linewidth', 1.2);
    plot(h,brydown,'--k','linewidth', 1.2);
    lgnd = legend(legendArr,...
        'Location','Best','Box','off');
        
    set(lgnd,'FontSize',14);
    xlim([0,1-J]);
    yticks(0.2:0.2:1.0);
    xticks(0:(1-J)/2:1-J);
    set(gca,'xtick',0:0.2:1.0);
    set(gca,'FontSize',14)
    xlabel('{\bf Flux factor} $h$','Interpreter','LaTeX','FontSize',19);
    ylabel('{\bf Eddington factor} $\chi$','Interpreter','LaTeX','FontSize',19);
    title(['$\mathcal{J}$ = ',num2str(J)],'Interpreter','LaTeX','FontSize',20);
    
    x0=10;
    y0=5;
    width=500;
    height=400;
    set(gcf,'units','points','position',[x0,y0,width,height])
    saveas(fig,['Closures0_', num2str(J*100),'.png']);
end

function [ value ] = Closure( J, h, theta_fun )

temp = 2.*(1-J).*(1-2.*J)/3 .* theta_fun(h/(1-J));
value = 1/3 + temp;

end
