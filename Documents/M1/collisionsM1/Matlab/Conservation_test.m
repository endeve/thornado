clear all

%%% including electron-type antineutrinos 
%%% including Pair opacity

readTestConstants();
readTables();

ReadMatterProfile;

global g_E_N g_W2_N g_W3_N
% generate neutrino energy grid
Numbins = 16;
NumPts = 2;
[g_E_N, g_W2_N, g_W3_N] = ComputePointsAndWeightsE(Numbins, NumPts);

global BoltzmannConstant

NumTests = size(MatterProfile,1);
% NumTests = 10;

Iters = zeros(NumTests,3);
Diffs =  zeros(NumTests,3);
DiffT = Diffs;
DiffY = Diffs;
DiffE = Diffs;

profile clear
profile on
for i = 1:NumTests

    status_text = ['Test ', num2str(i), ' out of ', num2str(NumTests)];
    disp(status_text);
    
    D = MatterProfile(i,2);
    T = MatterProfile(i,3);
    Y = MatterProfile(i,4);


    % dt is around 10^{-6} second
    dt = 1e-6;

    % compute equilibrium distribution 
    [Mnu, ~, ~] = ComputeNeutrinoChemicalPotentials(D, T, Y);
    [J0.Ne, J0.ANe] = FermiDirac( g_E_N, Mnu, BoltzmannConstant * T );
    
    % initial condition (perturbed equilibrium distribution)
%     J.Ne = J0.Ne - 1e-2*J0.Ne; 
%     J.ANe = J0.ANe - 1e-2*J0.ANe; 

    % initial condition (equilibrium distribution w/ perturbed (D,T,Y))
    J = J0;
    perturbation = 1e-1;
    D = D - perturbation * D;
    T = T - perturbation * T;
    Y = Y - perturbation * Y;

    % compute new equilibrium distribution 
    [Mnu, ~, ~] = ComputeNeutrinoChemicalPotentials(D, T, Y);
    [J0.Ne, J0.ANe] = FermiDirac( g_E_N, Mnu, BoltzmannConstant * T );
    
    % compute chi
    D_N = D * ones(size(g_E_N));
    T_N = T * ones(size(g_E_N));
    Y_N = Y * ones(size(g_E_N));

    [Chi.Ne, Chi.ANe] = ComputeAbEmOpacityOnEGrid_TABLE(g_E_N, D_N, T_N, Y_N);
    
    % compute internal energy
    [E, dEdT, dEdY] = ComputeSpecificInternalEnergy_TABLE(D, T, Y);

    % Pair kernel test
    [R_Pr_Ne, R_An_Ne, R_Pr_ANe, R_An_ANe] = ComputePairOpacityOnEGrid_TABLE(g_E_N, D, T, Y);
    
    % Newton's method
    % [J01, D1, T1, Y1, E1, iter1] = SolveMatterEquations_EmAb( J, dt * Chi, D, T, Y, E );
    % [J00, D0, T0, Y0, E0, iter0] = SolveMatterEquations_EmAb_orig( J, dt * Chi, D, T, Y, E );

    % Fixed-point
    % [J02, D2, T2, Y2, E2, iter2, Inneriter2] = SolveMatterEquations_EmAb_FP( J, dt * Chi, D, T, Y, E );
    % explicit NES
    % [J02, D2, T2, Y2, E2, iter2, Inneriter2] = SolveMatterEquations_EmAb_NES_FP( J, dt, Chi, R_in_NES, R_out_NES, D, T, Y, E );
    
    % nested FP
%     [J01, D1, T1, Y1, E1, iter1, Inneriter1] = SolveMatterEquations_EmAb_NES_FP( J, dt, Chi, D, T, Y, E );
%     [J01, D1, T1, Y1, E1, iter1AA, Inneriter1AA] = SolveMatterEquations_EmAb_NES_FP_AA( J, dt, Chi, D, T, Y, E );
    Inneriter1 = 0;
    Inneriter1AA = 0;
    iter1 = 0;
    iter1AA = 0;
    
    % coupled FP
%     [J02, D2, T2, Y2, E2, iter2backup] = SolveMatterEquations_EmAb_NES_FP_coupled_backup( J, dt, Chi, D, T, Y, E );
%     [J0_NES, D_NES, T_NES, Y_NES, E_NES, iter_NES] = SolveMatterEquations_EmAb_NES_FP_coupled( J, dt, Chi, D, T, Y, E );

%     [J0_NES, D_NES, T_NES, Y_NES, E_NES, iter_NES] = SolveMatterEquations_EmAb_NES_FP_coupledAA_noAN( J, dt, Chi, D, T, Y, E );
%     [J0_NESAN, D_NESAN, T_NESAN, Y_NESAN, E_NESAN, iter_NESAN] = SolveMatterEquations_EmAb_NES_FP_coupledAA( J, dt, Chi, D, T, Y, E );
    [J0_TP, D_TP, T_TP, Y_TP, E_TP, iter_TP] = SolveMatterEquations_Pair( J, dt, Chi, D, T, Y, E );

%     Iters(i,:) = [iter1 Inneriter1 iter1AA Inneriter1AA iter2backup iter2 iter2AA];
    Iters(i,:) = [iter_TP];
%     Iters(i,:) = [iter_NES iter_NESAN iter_TP];
%     Diffs(i,1) = norm(([T_NES Y_NES E_NES]-[T_NESAN Y_NESAN E_NESAN])./[T_NES Y_NES E_NES]);
%     Diffs(i,2) = norm(([T_NES Y_NES E_NES]-[T_TP Y_TP E_TP])./[T_NES Y_NES E_NES]);
%     Diffs(i,3) = norm(([T_TP Y_TP E_TP]-[T_NESAN Y_NESAN E_NESAN])./[T_NES Y_NES E_NES]);
%     DiffT(i,:) = abs([(T_NES - T_NESAN), (T_NES - T_TP), (T_NESAN - T_TP) ] ./ T_NES);
%     DiffY(i,:) = abs([(Y_NES - Y_NESAN), (Y_NES - Y_TP), (Y_NESAN - Y_TP) ] ./ Y_NES);
%     DiffE(i,:) = abs([(E_NES - E_NESAN), (E_NES - E_TP), (E_NESAN - E_TP) ] ./ E_NES);
    for k = 1:length(status_text)+1
        fprintf('\b');
    end
end

finish_text = ['All ', num2str(NumTests), ' tests completed.'];
disp(finish_text);

profile off
profile viewer