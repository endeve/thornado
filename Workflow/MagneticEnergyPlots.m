filesplusone = 102;

TimeVec           = zeros(filesplusone);
TotalMagEnergyVec = zeros(filesplusone);
TotalDivEnergyVec = zeros(filesplusone);

ProblemName = 'Advection2D';
Directory = '//wsl$/ubuntu/home/jbuffal/thornado_MHD_3D/SandBox/dgExperiments_MHD_Relativistic_IDEAL/Output/LoopAdvectionStatic/LoopAdvectionStatic_2ndOrder_n128x128x1_NoLimiting_DivClean_Damp0.00e0';
nNodes = 2;

[ Time, X1, X2, X3, CD, CS1, CS2, CS3, CE, CNe, CB1, CB2, CB3, CChi] = ...
    ReadMagnetofluidFields_Conserved( ProblemName, 0, Directory);
[ ~, ~, ~, ~, PD, PV1, PV2, PV3, PE, PNe, PB1, PB2, PB3, PChi] = ...
    ReadMagnetofluidFields_Primitive( ProblemName, 0, Directory);
[ ~, ~, ~, ~, AM_P, AM_T, AM_Ye, AM_S, AM_E, AM_Me, AM_Mp, AM_Mn, AM_Yp, AM_Yn, AM_Ya, AM_Yh, AM_Gm, AM_Cs ] = ...
    ReadMagnetofluidFields_Auxiliary( ProblemName, 0, Directory);

Dim = 1;
dX1 = X1(1 + nNodes) - X1(1);
dV = dX1;
if(size(X2,1) > 1)
    Dim = 2;
    dX2 = X2(1 + nNodes) - X2(1);
    dV = dX1 * dX2;
elseif(size(X3,1) > 1)
    Dim = 3;
    dX3 = X3(1 + nNodes) - X3(1);
    dV = dX1 * dX2 * dX3;
end

W = 1.0 ./ sqrt(1.0 - (PV1.^2 + PV2.^2 + PV3.^2));
b0 = W .* (CB1 .* PV1 + CB2 .* PV2 + CB3 .* PV3);
MagEnergy = 0.5 * ((CB1.^2 + CB2.^2 + CB3.^2) + b0.^2) ./ W.^2;
DivEnergy = 0.5 * CChi.^2;

MagEnergyAvg = GetCellAverages(Dim, nNodes, MagEnergy);
DivEnergyAvg = GetCellAverages(Dim, nNodes, DivEnergy);

TotalMagEnergy = sum(sum(sum(MagEnergyAvg * dV)));
TotalDivEnergy = sum(sum(sum(DivEnergyAvg * dV)));

TimeVec(1) = Time;
TotalMagEnergyVec(1) = TotalMagEnergy;
TotalDivEnergyVec(1) = TotalDivEnergy;

for k = 1:filesplusone - 1

    [ Time, X1, X2, X3, CD, CS1, CS2, CS3, CE, CNe, CB1, CB2, CB3, CChi] = ...
        ReadMagnetofluidFields_Conserved( ProblemName, k, Directory);
    [ ~, ~, ~, ~, PD, PV1, PV2, PV3, PE, PNe, PB1, PB2, PB3, PChi] = ...
        ReadMagnetofluidFields_Primitive( ProblemName, k, Directory);
    [ ~, ~, ~, ~, AM_P, AM_T, AM_Ye, AM_S, AM_E, AM_Me, AM_Mp, AM_Mn, AM_Yp, AM_Yn, AM_Ya, AM_Yh, AM_Gm, AM_Cs ] = ...
        ReadMagnetofluidFields_Auxiliary( ProblemName, k, Directory);     

    W = 1.0 ./ sqrt(1.0 - (PV1.^2 + PV2.^2 + PV3.^2));
    b0 = W .* (CB1 .* PV1 + CB2 .* PV2 + CB3 .* PV3);
    MagEnergy = 0.5 * ((CB1.^2 + CB2.^2 + CB3.^2) + b0.^2) ./ W.^2;
    DivEnergy = 0.5 * CChi.^2;

    MagEnergyAvg = GetCellAverages(Dim, nNodes, MagEnergy);
    DivEnergyAvg = GetCellAverages(Dim, nNodes, DivEnergy);

    TotalMagEnergy = sum(sum(sum(MagEnergyAvg * dV)));
    TotalDivEnergy = sum(sum(sum(DivEnergyAvg * dV)));
    
    TotalMagEnergyVec(k+1) = TotalMagEnergy;
    TotalDivEnergyVec(k+1) = TotalDivEnergy;
    TimeVec(k+1) = Time;

end
semilogy(TimeVec, TotalMagEnergyVec, 'linewidth', 2.5);
hold on
semilogy(TimeVec, TotalDivEnergyVec, 'linewidth', 2.5);