filesplusone = 102;
angle = pi/4.0;
VA = 0.38196601125010510;

TimeVec        = zeros(filesplusone);
MaxErrorVec    = zeros(filesplusone);
MaxRelErrorVec = zeros(filesplusone);

M(filesplusone) = struct('cdata', [], 'colormap', []);

Directory = '//wsl$/ubuntu/home/jbuffal/thornado_MHD_3D/SandBox/dgExperiments_MHD_Relativistic_IDEAL/Output/DivergenceCleaning/NonSmooth_1stOrder_n128x1x1_NoLimiting_DivClean_Damp0.00e00';
fname = 'C:\Users\jbuff\Desktop\Test.mp4';
delete  'C:\Users\jbuff\Desktop\Test.mp4';
Problem   = 'Cleaning1D';
NumFields = 1;

[ Time, X1, X2, X3, CD_M, CS1_M, CS2_M, CS3_M, CE_M, CNe_M, CB1_M, CB2_M, CB3_M, CChi_M] = ...
    ReadMagnetofluidFields_Conserved( Problem, 0, Directory);
[ ~, ~, ~, ~, PD_M, PV1_M, PV2_M, PV3_M, PE_M, PNe_M, PB1_M, PB2_M, PB3_M, PChi_M] = ...
    ReadMagnetofluidFields_Primitive( Problem, 0, Directory);
[ ~, ~, ~, ~, AP_M, AT_M, AYe_M, AS_M, AE_M, AMe_M, AMp_M, AMn_M, AYp_M, AYn_M, AYa_M, AYh_M, AGM_M, ACs_M ] = ...
    ReadMagnetofluidFields_Auxiliary( Problem, 0, Directory);        
[ ~, ~, ~, ~, DTCI_M, DShockX1_M, DShockX2_M, DShockX3_M, DTheta1_M, DTheta2_M, DTheta3_M, DMinE_M, DMaxE_M, DDiv_M ] = ...
    ReadMagnetofluidFields_Diagnostic( Problem, 0, Directory );

% AnalyticSolution = zeros(size(X1,1), size(X2,1));
InitialFields = {CB1_M'};
% for i = 1:size(X1,1)
%     for j = 1:size(X2,1)
%         AnalyticSolution(i,j) = VA * sin(angle)* cos(2 * pi * (X1(i) * cos(angle) + X2(j) * sin(angle)));
%     end
% end

%TimeVec(1) = Time;
%Error    = abs(Field - AnalyticSolution');
%RelError = Error ./ abs(AnalyticSolution');
%MaxErrorVec(1)    = max(max(Error));
%MaxRelErrorVec(1) = max(max(RelError));

h = figure;
h.Units = 'normalized';
h.Position = [0 0 0.5 1.0];

plot( X1, InitialFields{1}, 'linewidth', 2.5);

ax = gca;
ax.FontSize = 14;
ax.FontWeight = 'bold';
%legend('E1 - N1', 'E1 - N2', 'E1 - N3', 'E2 - N1', 'E2 - N2', 'E2 - N3');
title({['Cleaning 1D - Non-Smooth, 1st Order, n0128x1x1, No Limiting, Divergence Cleaning'] ['\kappa = 0.00e00']});
subtitle('t = ' + compose("%4.2e", Time));
xlabel('x', 'FontWeight', 'bold');
ylabel('$\nabla \cdot \vec{B}$', 'FontWeight', 'bold', 'Interpreter', 'latex');
xlim([-1 1]);
ylim([-2 2])
ax.NextPlot = 'replaceChildren';

v = VideoWriter(fname, 'MPEG-4');
v.Quality = 100;
v.FrameRate = filesplusone/20;
open(v);
M(1) = getframe(h);
writeVideo(v, M(1));

for k = 1:filesplusone - 1

[ Time, X1, X2, X3, CD_M, CS1_M, CS2_M, CS3_M, CE_M, CNe_M, CB1_M, CB2_M, CB3_M, CChi_M] = ...
    ReadMagnetofluidFields_Conserved( Problem, k, Directory);
[ ~, ~, ~, ~, PD_M, PV1_M, PV2_M, PV3_M, PE_M, PNe_M, PB1_M, PB2_M, PB3_M, PChi_M] = ...
    ReadMagnetofluidFields_Conserved( Problem, k, Directory);
[ ~, ~, ~, ~, AP_M, AT_M, AYe_M, AS_M, AE_M, AMe_M, AMp_M, AMn_M, AYp_M, AYn_M, AYa_M, AYh_M, AGM_M, ACs_M ] = ...
    ReadMagnetofluidFields_Auxiliary( Problem, k, Directory);    
[ ~, ~, ~, ~, DTCI_M, DShockX1_M, DShockX2_M, DShockX3_M, DTheta1_M, DTheta2_M, DTheta3_M, DMinE_M, DMaxE_M, DDiv_M ] = ...
    ReadMagnetofluidFields_Diagnostic( Problem, k, Directory );

    Fields = {CB1_M'};
    % for i = 1:size(X1,1)
    %   for j = 1:size(X2,1)
    %         AnalyticSolution(i,j) = VA * sin(angle)* cos(2 * pi * (X1(i) * cos(angle) + X2(j) * sin(angle) - VA * Time));
    %     end
    % end
    % 
    % Error    = abs(Field - AnalyticSolution');
    % RelError = Error ./ abs(AnalyticSolution');

    %TimeVec(k+1)        = Time;
    %MaxErrorVec(k+1)    = max(max(Error));
    %MaxRelErrorVec(k+1) = max(max(RelError));

    plot( X1, Fields{1}, 'linewidth', 2.5);
    %legend('E1 - N1', 'E1 - N2', 'E1 - N3', 'E2 - N1', 'E2 - N2', 'E2 - N3');
    
    subtitle('t = ' + compose("%4.2e", Time));
    drawnow
    M(k+1) = getframe(h);
    writeVideo(v, M(k+1));
end

close(v);