filesplusone = 102;

TimeVec           = zeros(filesplusone);
TotalMagEnergyVec    = zeros(filesplusone);
TotalDivEnergyVec = zeros(filesplusone);

ProblemName = 'Advection2D';
Directory = '//wsl$/ubuntu/home/jbuffal/thornado_MHD_3D/SandBox/dgExperiments_MHD_Relativistic_IDEAL/Output/LoopAdvectionStatic/LoopAdvectionStatic_1stOrder_n128x128x1_NoLimiting_DivClean_Damp1.00e0';
Order = 1;
Dim = 2;

[ Time, X1, X2, X3, CD, CS1, CS2, CS3, CE, CNe, CB1, CB2, CB3, CChi] = ...
    ReadMagnetofluidFields_Conserved( ProblemName, 0, Directory);
[ ~, ~, ~, ~, PD, PV1, PV2, PV3, PE, PNe, PB1, PB2, PB3, PChi] = ...
    ReadMagnetofluidFields_Primitive( ProblemName, 0, Directory);
[ ~, ~, ~, ~, AM_P, AM_T, AM_Ye, AM_S, AM_E, AM_Me, AM_Mp, AM_Mn, AM_Yp, AM_Yn, AM_Ya, AM_Yh, AM_Gm, AM_Cs ] = ...
    ReadMagnetofluidFields_Auxiliary( ProblemName, 0, Directory);    

dX1 = X1(1 + Order) - X1(1)
dX2 = X2(1 + Order) - X2(1)
W = 1.0 ./ sqrt(1.0 - (PV1.^2 + PV2.^2 + PV3.^2));
min(min(W.^2));
max(max(W.^2));
b0 = W .* (CB1 .* PV1 + CB2 .* PV2 + CB3 .* PV3);
MagEnergy = 0.5 * ((CB1.^2 + CB2.^2 + CB3.^2) + b0.^2) ./ W.^2;
DivEnergy = 0.5 * CChi.^2;

MagEnergyAvg = GetCellAverages(Order, MagEnergy);
DivEnergyAvg = GetCellAverages(Order, DivEnergy);

TotalMagEnergy = sum(sum(sum(MagEnergyAvg * dX1 * dX2)));
TotalDivEnergy = sum(sum(sum(DivEnergyAvg * dX1 * dX2)));

TimeVec(1) = Time;
TotalMagEnergyVec(1) = TotalMagEnergy;
TotalDivEnergyVec(1) = TotalDivEnergy;

for k = 1:filesplusone - 1

    k

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

    MagEnergyAvg = GetCellAverages(Order, MagEnergy);
    DivEnergyAvg = GetCellAverages(Order, DivEnergy);

    TotalMagEnergy = sum(sum(sum(MagEnergyAvg * dX1 * dX2)));
    TotalDivEnergy = sum(sum(sum(DivEnergyAvg * dX1 * dX2)));
    
    TotalMagEnergyVec(k+1) = TotalMagEnergy;
    TotalDivEnergyVec(k+1) = TotalDivEnergy;
    TimeVec(k+1) = Time;

end
plot(TimeVec, TotalMagEnergyVec, 'linewidth', 2.5);
hold on
plot(TimeVec, TotalDivEnergyVec, 'linewidth', 2.5);