clear all;
% close all;
plot_track = false;
saveData = true;
%%%%%%% 
% start from Gaussian distributions on both nu_e and nu_e_bar
% evolve the collision equation until equilibrium reached
%%%%%%%

readTestConstants();
readTables();

global g_E_N g_W2_N g_W3_N;
% generate neutrino energy grid
Numbins = 16;
NumPts = 2;
[g_E_N, g_W2_N, g_W3_N] = ComputePointsAndWeightsE(Numbins, NumPts);


global BoltzmannConstant AtomicMassUnit PlanckConstant SpeedOfLight Erg2MeV;

ReadMatterProfile;
NumTests = size(MatterProfile,1);
% NumTests = 10;

hc = PlanckConstant * SpeedOfLight;
c = SpeedOfLight;
Theta2_N = 4 * pi * g_W2_N;
Theta3_N = 4 * pi * g_W3_N;

Theta2_N = Theta2_N / (hc)^3;
Theta3_N = Theta3_N / (hc)^3;


% profile clear
% profile on

% for density = [1e8 1e10 1e12]
% for density = [1e8 1e10 1e12 1e14]

for density = [1e14]
    
test_idx = find(MatterProfile(:,2)>=density,1);
% test_idx = 1;
% test_idx = round(NumTests/2);

% test_idx = NumTests;
% test_idx = NumTests-68;



D0 = MatterProfile(test_idx,2);
T0 = MatterProfile(test_idx,3);
Y0 = MatterProfile(test_idx,4);

% % original Gaussian initial condition
Gau_mean = 2.0 * BoltzmannConstant * T0;
% % Gau_mean = 150;
Gau_std = 10;
Gau_dist = 0.99 * exp(-(g_E_N-Gau_mean).^2/2/Gau_std^2);
% Gau_dist = Gau_dist./max(Gau_dist)*0.99;
% % Gau_dist = Gau_dist./max(Gau_dist)*0.49;

% % % thin Gaussian
% Gau_mean = 10.0* BoltzmannConstant * T0;
% Gau_std = 1;
% Gau_dist = 1/Gau_std/2/sqrt(pi) * exp(-(g_E_N-Gau_mean).^2/2/Gau_std^2);
% Gau_dist = Gau_dist./max(Gau_dist)*0.99;



% compute chemical potential and derivatives
[Mnu, ~, ~] = ComputeNeutrinoChemicalPotentials(D0, T0, Y0);
% equilibrium distribution
[J00.Ne, J00.ANe] = FermiDirac( g_E_N, Mnu, BoltzmannConstant * T0 );


% for dt = [1e-5 1e-6 1e-7]
% for dt = [1e-5 1e-6]
for dt = [1e-6]
    
D = D0; T = T0; Y = Y0;
J.Ne = Gau_dist;
J.ANe = Gau_dist;



switch(dt)
    case 1e-4
        NumSteps = 20;
    
    case 1e-5
        NumSteps = 2000;
    case 1e-6
        NumSteps = 1000;
    case 1e-7
        NumSteps = 200000;
end
Iters = zeros(NumSteps,4);
Iters_in = zeros(NumSteps,3);
E_hist = zeros(NumSteps,2);
N_hist = zeros(NumSteps,2);
T_matter = zeros(NumSteps,1);
Y_hist = zeros(NumSteps,1);
Mnu_matter = zeros(NumSteps,1);
Lepton = zeros(NumSteps,1);
Energy = zeros(NumSteps,1);

% dt is around 10^{-6} second
% dt = 1e-7;
t = 0;
i = 0;
t_final = dt * NumSteps;

plot_idx = round(linspace(1,NumSteps,100));
plot_counter =1;

profile on
while(i<NumSteps)

    
    if ((plot_track) && (i==plot_idx(plot_counter)))
        plot_counter = plot_counter + 1;
        figure(1);
        subplot(3,1,1)
        plot(g_E_N, Gau_dist);
        legend('Initial distribution')
        title(['t = ', num2str(0, '%.1e'), ', D = ', num2str(D0, '%.3e'), ', T = ', num2str(T0, '%.3e'), ', Y = ', num2str(Y0, '%.3e')]);
        subplot(3,1,2)
        plot(g_E_N, J.Ne);
        %hold on
        legend('$J_{\nu_e}$', 'Interpreter', 'LaTex');
        title(['t = ', num2str(t, '%.1e'), ', D = ', num2str(D, '%.3e'), ', T = ', num2str(T, '%.3e'), ', Y = ', num2str(Y, '%.3e')]);
        subplot(3,1,3)
        plot(g_E_N, J.ANe);
        %hold on
        legend('$\bar{J}_{\nu_e}$', 'Interpreter', 'LaTex');
        title(['t = ', num2str(t, '%.1e'), ', D = ', num2str(D, '%.3e'), ', T = ', num2str(T, '%.3e'), ', Y = ', num2str(Y, '%.3e')]);
        pause(0.1);
    end
    i = i + 1;
   
    status_text = ['Time_step = ', num2str(i), ' out of ', num2str(NumSteps)];
    disp(status_text);
    
    % compute chi
    D_N = D * ones(size(g_E_N));
    T_N = T * ones(size(g_E_N));
    Y_N = Y * ones(size(g_E_N));

    [Chi.Ne, Chi.ANe] = ComputeAbEmOpacityOnEGrid_TABLE(g_E_N, D_N, T_N, Y_N);
    
    % compute internal energy
    [E, dEdT, dEdY] = ComputeSpecificInternalEnergy_TABLE(D, T, Y);

    % solve collision equation
     [~, ~, D_Coupled, T_Coupled, Y_Coupled, E_Coupled, iter_TP] = SolveMatterEquations_Pair( J, dt, Chi, D, T, Y, E );
    iter_in_TP = 0;
    [~, ~, ~, ~, ~, ~, iter_TPPI, iter_in_PI] = SolveMatterEquations_Pair_Nested( J, dt, Chi, D, T, Y, E );
    [~, ~, ~, ~, ~, ~, iter_TPAA, iter_in_AA] = SolveMatterEquations_Pair_NestedAA( J, dt, Chi, D, T, Y, E );
    [J_TP, J0_TP, D_TP, T_TP, Y_TP, E_TP, iter_TPNewton, iter_in_Newton] = SolveMatterEquations_Pair_NestedNewton( J, dt, Chi, D, T, Y, E );

    diff = abs(T_TP - T_Coupled)/abs(T_Coupled) + abs(Y_TP - Y_Coupled)/abs(Y_Coupled) + abs(E_TP - E_Coupled)/abs(E_Coupled);
    
    % update variables
    D = D_TP; T = T_TP; Y = Y_TP;
    
    J = J_TP;

    J0 = J0_TP;
    

    % store data
    N_B = D / AtomicMassUnit;

    Lepton(i) = sum(Theta2_N.*(J.Ne - J.ANe)) + N_B * Y ; %scaled
    Energy(i) = sum(Theta3_N.*(J.Ne + J.ANe)) + D * E_TP * Erg2MeV ; %scaled
    N_hist(i,1) = sum(g_W2_N.*J.Ne);
    E_hist(i,1) = sum(g_W3_N.*J.Ne);
    N_hist(i,2) = sum(g_W2_N.*J.ANe);
    E_hist(i,2) = sum(g_W3_N.*J.ANe);
    T_matter(i) = T;
    Y_hist(i) = Y;
    [Mnu, ~, ~] = ComputeNeutrinoChemicalPotentials(D, T, Y);
    Mnu_matter(i) = Mnu;
    Iters(i,:) = [iter_TPPI iter_TPAA iter_TPNewton iter_TP];
    Iters_in(i,:) = [iter_in_PI iter_in_AA iter_in_Newton];

    
    % evolve time
    t = t + dt;
    
    for k = 1:length(status_text)+1
        fprintf('\b');
    end
    
end

finish_text = ['All ', num2str(NumSteps), ' steps completed.'];
disp(finish_text);

profile off
profile viewer

% figure;
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
% 
% 
% figure;
% subplot(4,1,1)
% plot(g_E_N, J00.Ne);
% legend('$J_{0,\nu_e}$', 'Interpreter', 'LaTex');
% title(['t = ', num2str(0, '%.1e'), ', D = ', num2str(D0, '%.3e'), ', T = ', num2str(T0, '%.3e'), ', Y = ', num2str(Y0, '%.3e')]);
% subplot(4,1,2)
% plot(g_E_N, J00.ANe);
% legend('$\bar{J}_{0,\nu_e}$', 'Interpreter', 'LaTex');
% title(['t = ', num2str(0, '%.1e'), ', D = ', num2str(D0, '%.3e'), ', T = ', num2str(T0, '%.3e'), ', Y = ', num2str(Y0, '%.3e')]);
% subplot(4,1,3)
% plot(g_E_N, J0.Ne);
% legend('$J_{0,\nu_e}$', 'Interpreter', 'LaTex');
% title(['t = ', num2str(t, '%.1e'), ', D = ', num2str(D, '%.3e'), ', T = ', num2str(T, '%.3e'), ', Y = ', num2str(Y, '%.3e')]);
% subplot(4,1,4)
% plot(g_E_N, J0.ANe);
% legend('$\bar{J}_{0,\nu_e}$', 'Interpreter', 'LaTex');
% title(['t = ', num2str(t, '%.1e'), ', D = ', num2str(D, '%.3e'), ', T = ', num2str(T, '%.3e'), ', Y = ', num2str(Y, '%.3e')]);
% 
% figure;
% subplot(1,2,1)
% semilogx(Iters(:,1));
% hold on
% semilogx(Iters(:,2));
% semilogx(Iters(:,3));
% semilogx(Iters(:,4));
% title(['Iteration count, ''Density = ', num2str(density,'%.0e')])
% legend('PI','AA','Newton', 'Coupled')
% % semilogx(ones(1,NumSteps),':');
% 
% subplot(1,2,2)
% semilogx(Iters_in(:,1));
% hold on
% semilogx(Iters_in(:,2));
% semilogx(Iters_in(:,3));
% title(['[Inner iteration count, ''Density = ', num2str(density,'%.0e')])
% legend('PI','AA','Newton')
% % semilogx(ones(1,NumSteps),':');


if(saveData)
    save(['Data_' num2str(test_idx) '_dt_' num2str(dt, '%.0e') '_step_' num2str(NumSteps,'%.0e') '_thin_Gaussian_new.mat'],...
        'D0', 'T0', 'Y0', 'Gau_dist', 'J00', 'D', 'T', 'Y', 'J', 'J0', 'Iters', 'Iters_in', 'g_E_N', ...
        'E_hist', 'N_hist', 'T_matter', 'Y_hist', 'Mnu_matter', 'Lepton', 'Energy');
end
end
end
