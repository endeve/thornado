filesplusone = 102;
angle = pi/4.0;
VA = 0.38196601125010510;

Directory = ['//wsl$/ubuntu/home/jbuffal/thornado_MHD_3D/SandBox/dgExperiments_MHD_Relativistic_IDEAL/Output/OrszagTang2D/OT2D_1stOrder_n64x64x1_NoLimiting_CleaningOn_Damp0.00e00'];
SaveFile  =  ['C:\Users\jbuff\Desktop\Research\MHD\Numerical Results - New\OrszagTang\Movies\OT_1stOrder_n64x64x1_NoLimiting_CleaningOn_Damp0.00e00_PD.mp4'];
Problem   = 'OrszagTang2D'

TimeVec        = zeros(filesplusone);
MaxErrorVec    = zeros(filesplusone);
MaxRelErrorVec = zeros(filesplusone);

M(filesplusone) = struct('cdata', [], 'colormap', []);

[ Time, X1, X2, X3, CD_M, CS1_M, CS2_M, CS3_M, CE_M, CNe_M, CB1_M, CB2_M, CB3_M, CChi_M] = ...
    ReadMagnetofluidFields_Conserved( Problem, 0, Directory);
[ ~, ~, ~, ~, PD_M, PV1_M, PV2_M, PV3_M, PE_M, PNe_M, PB1_M, PB2_M, PB3_M, PChi_M] = ...
    ReadMagnetofluidFields_Primitive( Problem, 0, Directory);
[ ~, ~, ~, ~, AP_M, AT_M, AYe_M, AS_M, AE_M, AMe_M, AMp_M, AMn_M, AYp_M, AYn_M, AYa_M, AYh_M, AGM_M, ACs_M ] = ...
    ReadMagnetofluidFields_Auxiliary( Problem, 0, Directory);    
[ ~, ~, ~, ~, DTCI_M, DShockX1_M, DShockX2_M, DShockX3_M, DTheta1_M, DTheta2_M, DTheta3_M, DMinE_M, DMaxE_M, DDiv_M ] = ...
    ReadMagnetofluidFields_Diagnostic( Problem, 0, Directory );

InitialField = PD_M';
%InitialField = log10( abs( ( PV3_M' - ( 1.0 / 24.0 ) ) ) / ( 1.0 / 24.0 ) );

h = figure;
h.Units = 'normalized';
h.Position = [0 0 0.5 1.0];

[X, Y] = meshgrid(X1, X2);
imagesc( X1, X2, InitialField);
ax = gca;
ax.YDir = 'normal';
ax.FontSize = 14;
ax.FontWeight = 'bold';
axis equal;
title({['Orszag-Tang, 1st Order, n64x64, No Limiting, Divergence Cleaning'] ['\kappa = 0.00e00']});
subtitle('t = ' + compose("%4.2e", Time));
xlabel('x', 'FontWeight', 'bold');
ylabel('y', 'FontWeight', 'bold');
xlim([0.0 1.0])
ylim([0.0 1.0])
colormap(copper);
cbar = colorbar;
cbar.FontSize = 14;
%clim([min(min(InitialField)) max(max(InitialField))])
%clim([cos(angle)-sin(angle) cos(angle)+sin(angle)])
%clim([-VA * sin(angle) VA * sin(angle)]);
clim([0 0.5])
cbar.Label.String = '$\rho$ (Cell Averages)';
cbar.Label.FontWeight = 'bold';
%cbar.Label.Rotation = 0;
cbar.Label.Interpreter = 'latex';

ax.NextPlot = 'replaceChildren';

delete 'C:\Users\jbuff\Desktop\Research\MHD\Numerical Results - New\OrszagTang\Movies\OT_1stOrder_n64x64x1_NoLimiting_CleaningOn_Damp0.00e00_PD.mp4';

v = VideoWriter(SaveFile, 'MPEG-4');
v.Quality = 100;
v.FrameRate = filesplusone/10;
open(v);
M(1) = getframe(h);
writeVideo(v, M(1));

for k = 1:filesplusone - 1

[ Time, X1, X2, X3, CD_M, CS1_M, CS2_M, CS3_M, CE_M, CNe_M, CB1_M, CB2_M, CB3_M, CChi_M] = ...
    ReadMagnetofluidFields_Conserved( Problem, k, Directory);
[ ~, ~, ~, ~, PD_M, PV1_M, PV2_M, PV3_M, PE_M, PNe_M, PB1_M, PB2_M, PB3_M, PChi_M] = ...
    ReadMagnetofluidFields_Primitive( Problem, k, Directory);
[ ~, ~, ~, ~, AP_M, AT_M, AYe_M, AS_M, AE_M, AMe_M, AMp_M, AMn_M, AYp_M, AYn_M, AYa_M, AYh_M, AGM_M, ACs_M ] = ...
    ReadMagnetofluidFields_Auxiliary( Problem, k, Directory);    
[ ~, ~, ~, ~, DTCI_M, DShockX1_M, DShockX2_M, DShockX3_M, DTheta1_M, DTheta2_M, DTheta3_M, DMinE_M, DMaxE_M, DDiv_M ] = ...
    ReadMagnetofluidFields_Diagnostic( Problem, k, Directory );

    Field = PD_M';
    %Field = log10( abs( ( PV3_M' - ( 1.0 / 24.0 ) ) ) / ( 1.0 / 24.0 ) );

    imagesc( X1, X2, Field);
    subtitle('t = ' + compose("%4.2e", Time));
    drawnow
    M(k+1) = getframe(h);
    writeVideo(v, M(k+1));
end

close(v);