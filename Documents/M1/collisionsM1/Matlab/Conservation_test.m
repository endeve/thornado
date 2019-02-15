clear all

% D = 7.833e10;   %  1.6605e+03 to 3.1641e+15
% T = 2.4447e11;  %  1.1605e+09 to 1.8392e+12
% Y = 1.421e-1;   %  0.01       to 0.6
D = 4.146E14;
T = 1.260E11;
Y = 2.518E-01;
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
Numbins = 20;
NumPts = 2;
[g_E_N, g_W2_N, g_W3_N] = ComputePointsAndWeightsE(Numbins, NumPts);


D_N = D * ones(size(g_E_N));
T_N = T * ones(size(g_E_N));
Y_N = Y * ones(size(g_E_N));

[Chi] = ComputeAbEmOpacityOnEGrid_TABLE(g_E_N, D_N, T_N, Y_N);



dt = 1e-6;
global BoltzmannConstant
J0 = FermiDirac( g_E_N, Mnu, BoltzmannConstant * T );

J = zeros(size(g_E_N));
J(1:end-1) = J0(1:end-1);
% J(1) = 1;

global g_Iterations_Min g_Iterations_Max g_Iterations_Ave;
g_Iterations_Min = 0;
g_Iterations_Max = Inf;
g_Iterations_Ave = 0;

SolveMatterEquations_EmAb( J, dt * Chi, D, T, Y, E );

% dt is around 10^{-6} second
