clear all

% D = 7.833e10;   %  1.6605e+03 to 3.1641e+15
% T = 2.4447e10;  %  1.1605e+09 to 1.8392e+12
% Y = 1.421e-1;   %  0.01       to 0.6
% profile = [  
%   1.179E+05   4.117E+14   1.261E+11   2.520E-01  
%   1.658E+06   2.469E+13   2.517E+11   1.849E-01
%   2.120E+07   1.107E+08   7.010E+09   4.979E-01  
%   8.545E+06   1.836E+10   4.068E+10   2.530E-01  
%   1.566E+06   2.837E+13   2.641E+11   1.979E-01  
%   2.901E+04   4.146E+14   1.260E+11   2.518E-01  
% ];

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

Iters = zeros(NumTests,3);
Diffs =  zeros(NumTests,1);

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
    J0 = FermiDirac( g_E_N, Mnu, BoltzmannConstant * T );
    
    % initial condition (perturbed equilibrium distribution)
%     J = J0 - 1e-2*J0; 

    % initial condition (equilibrium distribution w/ perturbed (D,T,Y))
    J = J0; 

    D = D - 1e-3 * D;
    T = T - 1e-3 * T;
    Y = Y - 1e-3 * Y;

    % compute new equilibrium distribution 
    [Mnu, ~, ~] = ComputeNeutrinoChemicalPotentials(D, T, Y);
    J0 = FermiDirac( g_E_N, Mnu, BoltzmannConstant * T );
    
    % compute chi
    D_N = D * ones(size(g_E_N));
    T_N = T * ones(size(g_E_N));
    Y_N = Y * ones(size(g_E_N));

    [Chi] = ComputeAbEmOpacityOnEGrid_TABLE(g_E_N, D_N, T_N, Y_N);
    
    % compute internal energy
    [E, dEdT, dEdY] = ComputeSpecificInternalEnergy_TABLE(D, T, Y);


    % Newton's method
    % [J01, D1, T1, Y1, E1, iter1] = SolveMatterEquations_EmAb( J, dt * Chi, D, T, Y, E );
    % [J00, D0, T0, Y0, E0, iter0] = SolveMatterEquations_EmAb_orig( J, dt * Chi, D, T, Y, E );

    % Fixed-point
    % [J02, D2, T2, Y2, E2, iter2, Inneriter2] = SolveMatterEquations_EmAb_FP( J, dt * Chi, D, T, Y, E );
    % explicit NES
    % [J02, D2, T2, Y2, E2, iter2, Inneriter2] = SolveMatterEquations_EmAb_NES_FP( J, dt, Chi, R_in_NES, R_out_NES, D, T, Y, E );
    
    % nested FP
    [J01, D1, T1, Y1, E1, iter1, Inneriter1] = SolveMatterEquations_EmAb_NES_FP( J, dt, Chi, D, T, Y, E );
    
    % coupled FP
    [J02, D2, T2, Y2, E2, iter2] = SolveMatterEquations_EmAb_NES_FP_coupled( J, dt, Chi, D, T, Y, E );
    
    Iters(i,:) = [iter1 Inneriter1 iter2];
    Diffs(i) = norm(([T1 Y1 E1]-[T2 Y2 E2])./[T1 Y1 E1]);
    
    for k = 1:length(status_text)+1
        fprintf('\b');
    end
end

finish_text = ['All ', num2str(NumTests), ' tests completed.'];
disp(finish_text);
