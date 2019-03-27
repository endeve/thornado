clear all

% D = 7.833e10;   %  1.6605e+03 to 3.1641e+15
% T = 2.4447e10;  %  1.1605e+09 to 1.8392e+12
% Y = 1.421e-1;   %  0.01       to 0.6

profile = [  
  1.179E+05   4.117E+14   1.261E+11   2.520E-01  
];

D = profile(2);
T = profile(3);
Y = profile(4);

% D = flip(logspace(3.5, 15, 100)');
% T = (logspace(9.5, 12, 100)');
% Y = (logspace(-1.9, -0.9, 100)');

% readConstants();  
readTestConstants();
readTables();
[Mnu, ~, ~] = ComputeNeutrinoChemicalPotentials(D, T, Y);
% readtestTables();
[E, dEdT, dEdY] = ComputeSpecificInternalEnergy_TABLE(D, T, Y);


% T_recons = ComputeTemperatureFromSpecificInternalEnergy_TABLE(D, E, Y);
% norm((T-T_recons)./T)
    
% readTestConstants(); 

global g_E_N g_W2_N g_W3_N
% generate neutrino energy grid
Numbins = 16;
NumPts = 2;
[g_E_N, g_W2_N, g_W3_N] = ComputePointsAndWeightsE(Numbins, NumPts);


D_N = D * ones(size(g_E_N));
T_N = T * ones(size(g_E_N));
Y_N = Y * ones(size(g_E_N));

[Chi] = ComputeAbEmOpacityOnEGrid_TABLE(g_E_N, D_N, T_N, Y_N);

% explicit NES test
[R_in_NES, R_out_NES] = ComputeNesScatteringOpacityOnEGrid_TABLE(g_E_N, D, T, Y);

% dt is around 10^{-6} second
dt = 1e-6;

% compute equilibrium distribution 
global BoltzmannConstant
J0 = FermiDirac( g_E_N, Mnu, BoltzmannConstant * T );

J = zeros(size(g_E_N));
% J(1:end-1) = J0(1:end-1); % initial condition
% J(1) = 1;
J = J0 - 1e-2*J0; % initial condition

global g_Iterations_Min g_Iterations_Max g_Iterations_Ave;
g_Iterations_Min = 0;
g_Iterations_Max = Inf;
g_Iterations_Ave = 0;

% Newton's method
% [J01, D1, T1, Y1, E1, iter1] = SolveMatterEquations_EmAb( J, dt * Chi, D, T, Y, E );
% [J00, D0, T0, Y0, E0, iter0] = SolveMatterEquations_EmAb_orig( J, dt * Chi, D, T, Y, E );

% Fixed-point
% [J02, D2, T2, Y2, E2, iter2, Inneriter2] = SolveMatterEquations_EmAb_FP( J, dt * Chi, D, T, Y, E );
[J02, D2, T2, Y2, E2, iter2, Inneriter2] = SolveMatterEquations_EmAb_NES_FP( J, dt, Chi, R_in_NES, R_out_NES, D, T, Y, E );

