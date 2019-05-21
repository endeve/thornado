readTestConstants();
readTables();

global BoltzmannConstant SpeedOfLight;
global g_E_N g_W2_N g_W3_N;
Numbins = 16;
NumPts = 2;
[g_E_N, g_W2_N, g_W3_N] = ComputePointsAndWeightsE(Numbins, NumPts);
c = SpeedOfLight;

ReadMatterProfile;


Fig2 = figure;

% for density = [1e8 1e10 1e12 1e14]
for density = [1e12]

test_idx = find(MatterProfile(:,2)>=density,1);
    
for dt = [1e-7 1e-6 1e-5]
% for dt = [1e-1]
    
    switch(dt)
        case 1e-1
            NumSteps = 2000;
        case 1e-5
            NumSteps = 100;
        case 1e-6
            NumSteps = 1000;
        case 1e-7 
            NumSteps = 10000;
    end
% dt = 1e-5; NumSteps = 1000;
% dt = 1e-6; NumSteps = 10000;
% dt = 1e-7; NumSteps = 100000;

time =  dt*(1:NumSteps);


% load(['Data_' num2str(test_idx) '_dt_' num2str(dt, '%.0e') '_step_' num2str(NumSteps,'%.0e') '_.mat']);
% load(['Data_' num2str(test_idx) '_dt_' num2str(dt, '%.0e') '_step_' num2str(NumSteps,'%.0e') '_thin_Gaussian.mat']);
load(['Data_' num2str(test_idx) '_dt_' num2str(dt, '%.0e') '_step_' num2str(NumSteps,'%.0e') '_Conservation.mat']);

[Mnu, ~, ~] = ComputeNeutrinoChemicalPotentials(D, T, Y);
kT = BoltzmannConstant * T;
% 
% exponent = min( max( g_E_N - Mu / kT, -100 ), 100 );
% F = 1.0 ./ ( exp( exponent ) + 1 );
% F2 = sum(g_W2_M .* F) * kT^3;
% F3 = sum(g_W3_M .* F) * kT^4;

% x = linspace(0,8,100)';
% m = 1:100;
% F2 = -sum((-exp(x).^m)./(m.^3),2);
% F3 = -sum((-exp(x).^m)./(m.^4),2);


x = linspace(-10,8,1000);
% x = Mnu/kT;
exponent = min( max( repmat(g_E_N,1,length(x)) - repmat(x,length(g_E_N),1), -100 ), 100 );
F = 1.0 ./ ( exp( exponent ) + 1 );
F2 = sum(g_W2_N .* F);
F3 = sum(g_W3_N .* F);

a = E_hist(:,1); b = N_hist(:,1);

LHS = a./b.^(4/3);
RHS = F3./F2.^(4/3);

eta = zeros(length(LHS),1);
y2 = zeros(length(LHS),1);
y3 = zeros(length(LHS),1);
for i = 1:length(LHS)
[val, idx] = min(abs(LHS(i) - RHS));

eta(i) = x(idx);
y2(i) = (N_hist(end,1)/F2(idx))^(1/3);
y3(i) = (E_hist(end,1)/F3(idx))^(1/4);
end
Mnu_approx = eta.*y2;
T_approx = y2/BoltzmannConstant;

x = linspace(-10,8,1000);
exponent = min( max( repmat(g_E_N,1,length(x)) + repmat(x,length(g_E_N),1), -100 ), 100 );
F_bar = 1.0 ./ ( exp( exponent ) + 1 );
F2_bar = sum(g_W2_N .* F_bar);
F3_bar = sum(g_W3_N .* F_bar);

a = E_hist(:,2); b = N_hist(:,2);

LHS = a./b.^(4/3);
RHS = F3_bar./F2_bar.^(4/3);

eta_bar = zeros(length(LHS),1);
y2_bar = zeros(length(LHS),1);
y3_bar = zeros(length(LHS),1);
for i = 1:length(LHS)
[val, idx] = min(abs(LHS(i) - RHS));

eta_bar(i) = x(idx);
y2_bar(i) = (N_hist(end,2)/F2_bar(idx))^(1/3);
y3_bar(i) = (E_hist(end,2)/F3_bar(idx))^(1/4);
end
Mnubar_approx = eta_bar.*y2_bar;
Tbar_approx = y2_bar/BoltzmannConstant;




% Fig = figure;
% sgtitle(['$\rho$ = ', num2str(density,'%.0e'), ', $\quad\Delta t$ = ', num2str(dt,'%.0e')],'interpreter','latex')
% subplot(4,1,1)
% semilogx(time, Mnu_approx, 'linewidth', 2)
% hold on
% semilogx(time, Mnubar_approx, 'linewidth', 2)
% semilogx(time, Mnu_matter,'k', 'linewidth', 2)
% legend('Neutrino','Antineutrino','Matter')
% title('Chemical potential $\mu_\nu$','interpreter','latex')
% subplot(4,1,2)
% semilogx(time, T_approx, 'linewidth', 2)
% hold on
% semilogx(time, Tbar_approx, 'linewidth', 2)
% semilogx(time, T_matter,'k', 'linewidth', 2)
% legend('Neutrino','Antineutrino','Matter','location','southeast')
% title('Temprature $T$','interpreter','latex')
% 
% subplot(4,1,3)
% semilogx(time,Iters(:,4), 'linewidth', 2);
% hold on
% semilogx(time,Iters(:,2), 'linewidth', 2);
% semilogx(time,Iters(:,3), 'linewidth', 2);
% ylim([0 15])
% title('Outer iteration','interpreter','latex')
% legend('Coupled - AA','Nested - AA/AA','Nested - AA/Newton')
% % semilogx(ones(1,NumSteps),':');
% 
% co = get(gca,'colororder');
% set(gca,'colororder',co(2:end,:));
% subplot(4,1,4)
% semilogx(time,Iters_in(:,2), 'linewidth', 2, 'color', [0.8500 0.3250 0.0980]);
% hold on
% semilogx(time,Iters_in(:,3), 'linewidth', 2, 'color', [0.9290 0.6940 0.1250]);
% ylim([0 8])
% title('Inner iteration','interpreter','latex')
% legend('AA','Newton')

% % plot conserved quantities
figure(Fig2)
sgtitle('Relative deviation of conserved quantities','interpreter','latex')
subplot(2,1,1)
hold on
semilogx(time, (Lepton - Lepton(1))./abs(Lepton(1)), 'linewidth', 2)
title('Lepton number','interpreter','latex')
legend('$\Delta t = 10^{-7}$', '$\Delta t = 10^{-6}$', '$\Delta t = 10^{-5}$','interpreter','latex')
xlabel('time','interpreter','latex')
ylim([-7e-7 7e-7])
subplot(2,1,2)
hold on
semilogx(time, (Energy - Energy(1))./abs(Energy(1)), 'linewidth', 2)
title('Energy (MeV)','interpreter','latex')
legend('$\Delta t = 10^{-7}$', '$\Delta t = 10^{-6}$', '$\Delta t = 10^{-5}$','interpreter','latex')
xlabel('time','interpreter','latex')
ylim([-5e-8 5e-8])
 
% Fig3 = figure;
% subplot(3,1,1)
% plot(g_E_N, Gau_dist);
% legend('Initial distribution')
% title(['t = ', num2str(0, '%.1e'), ', D = ', num2str(D0, '%.3e'), ', T = ', num2str(T0, '%.3e'), ', Y = ', num2str(Y0, '%.3e')]);
% subplot(3,1,2)
% plot(g_E_N, J.Ne);
% legend('$J_{\nu_e}$', 'Interpreter', 'LaTex');
% title(['t = ', num2str(t, '%.1e'), ', D = ', num2str(D, '%.3e'), ', T = ', num2str(T, '%.3e'), ', Y = ', num2str(Y, '%.3e')]);
% subplot(3,1,3)
% plot(g_E_N, J.ANe);
% legend('$\bar{J}_{\nu_e}$', 'Interpreter', 'LaTex');
% title(['t = ', num2str(t, '%.1e'), ', D = ', num2str(D, '%.3e'), ', T = ', num2str(T, '%.3e'), ', Y = ', num2str(Y, '%.3e')]);
% 

% filename = ['Rho_' num2str(density, '%.0e') '_dt_' num2str(dt, '%.0e') ];
% set(Fig, 'Position',[200 100 1060 882], 'PaperPositionMode','auto');
% savefig(Fig,['Figures/', filename, '.fig']);
% print(Fig,['Figures/', filename],'-depsc');
end

end

filename = ['Rho_' num2str(density, '%.0e') '_Conservation' ];
set(Fig2, 'Position',[200 100 1060 882], 'PaperPositionMode','auto');
savefig(Fig2,['Figures/', filename, '.fig']);
print(Fig2,['Figures/', filename],'-depsc');

% 
% figure;
% subplot(2,1,1)
% semilogx(time, Mnu_approx, 'linewidth', 2)
% hold on
% semilogx(time, Mnubar_approx, 'linewidth', 2)
% semilogx(time, Mnu_matter,'k', 'linewidth', 2)
% 
% subplot(2,1,2)
% semilogx(time, T_approx, 'linewidth', 2)
% hold on
% semilogx(time, Tbar_approx, 'linewidth', 2)
% semilogx(time, T_matter,'k', 'linewidth', 2)
% 
% 
% figure;
% sgtitle(['Iteration counts, ' 'Density = ', num2str(density,'%.0e')])
% 
% subplot(1,2,1)
% semilogx(time,Iters(:,4), 'linewidth', 2);
% hold on
% semilogx(time,Iters(:,2), 'linewidth', 2);
% semilogx(time,Iters(:,3), 'linewidth', 2);
% title('Outer iteration')
% legend('Coupled','Nested - AA/AA',' Nested - AA/Newton')
% % semilogx(ones(1,NumSteps),':');
% 
% subplot(1,2,2)
% semilogx(time,Iters_in(:,2), 'linewidth', 2);
% hold on
% semilogx(time,Iters_in(:,3), 'linewidth', 2);
% title('Inner iteration')
% legend('AA','Newton')



% figure;
% subplot(2,1,1)
% plot(eta)
% hold on
% plot(eta_bar)
% plot(Mnu/kT*ones(length(eta),1))
% 
% subplot(2,1,2)
% plot(y2)
% hold on
% plot(y2_bar)
% plot(kT*ones(length(eta),1))

