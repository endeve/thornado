function [ TimeVec, TotalMagEnergyVec, TotalThermEnergyVec, TotalDivEnergyVec, TotalDivVec ]...
    = TotalPlots( filesplusone, ProblemName, Directory, nNodes )

TimeVec             = zeros([filesplusone 1]);
TotalMagEnergy      = 0.0;
TotalDivEnergy      = 0.0;
TotalThermEnergy    = 0.0;
TotalKinEnergy      = 0.0;
TotalDiv            = 0.0;
TotalAllEnergy      = 0.0;
TotalMagEnergyVec   = zeros([filesplusone 1]);
TotalDivVec         = zeros([filesplusone 1]);
TotalDivEnergyVec   = zeros([filesplusone 1]);
TotalThermEnergyVec = zeros([filesplusone 1]);
TotalKinEnergyVec   = zeros([filesplusone 1]);
TotalAllEnergyVec   = zeros([filesplusone 1]);

[ Time, X1, X2, X3, CD, CS1, CS2, CS3, CE, CNe, CB1, CB2, CB3, CChi] = ...
    ReadMagnetofluidFields_Conserved( ProblemName, 0, Directory);
[ ~, ~, ~, ~, PD, PV1, PV2, PV3, PE, PNe, PB1, PB2, PB3, PChi] = ...
    ReadMagnetofluidFields_Primitive( ProblemName, 0, Directory);
[ ~, ~, ~, ~, AP, AT, AYe, AS, AE, AMe, AMp, AMn, AYp, AYn, AYa, AYh, AGm, ACs ] = ...
    ReadMagnetofluidFields_Auxiliary( ProblemName, 0, Directory);
[ ~, ~, ~, ~, DTCI, DShockX1, DShockX2, DShockX3, DTheta1, DTheta2, DTheta3, DMinE, DMaxE, DDiv ]...
    = ReadMagnetofluidFields_Diagnostic( ProblemName, 0, Directory );

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
bu0 = W .* (CB1 .* PV1 + CB2 .* PV2 + CB3 .* PV3);
bd0 = -bu0;
PBSqr = bu0 .* bd0 + PB1.^2 + PB2.^2 + PB3.^2;

MagEnergy = (W.^2 - 1.0/2.0) .* PBSqr - (bu0).^2;
Div    = DDiv;
DivEnergy   = 0.5 * CChi.^2;
ThermEnergy = PE;
KinEnergy   = PD .* W.^2 + AGm .* PE .* W.^2;
AllEnergy   = MagEnergy + DivEnergy + ThermEnergy;

MagEnergyAvg     = GetCellAverages(Dim, nNodes, MagEnergy);
DivAvg           = abs(GetCellAverages(Dim, nNodes, Div));
DivEnergyAvg     = GetCellAverages(Dim, nNodes, DivEnergy);
ThermEnergyAvg   = GetCellAverages(Dim, nNodes, ThermEnergy);
KinEnergyAvg     = GetCellAverages(Dim, nNodes, KinEnergy);
AllEnergyAvg     = GetCellAverages(Dim, nNodes, AllEnergy);

TotalMagEnergy   = sum(MagEnergyAvg   * dV, "all");
TotalDiv         = sum(DivAvg         * dV, "all");
TotalDivEnergy   = sum(DivEnergyAvg   * dV, "all");
TotalThermEnergy = sum(ThermEnergyAvg * dV, "all");
TotalKinEnergy   = sum(KinEnergyAvg   * dV, "all");
TotalAllEnergy   = sum(AllEnergyAvg   * dV, "all");

TotalMagEnergyVec(1)   = TotalMagEnergy;
TotalDivVec(1)         = TotalDiv;
TotalDivEnergyVec(1)   = TotalDivEnergy;
TotalThermEnergyVec(1) = TotalThermEnergy;
TotalKinEnergyVec(1)   = TotalKinEnergy;
TotalAllEnergyVec(1)   = TotalAllEnergy;

TimeVec(1) = Time;

for k = 1:filesplusone - 1

    TotalMagEnergy    = 0.0;
    TotalDiv          = 0.0;
    TotalDivEnergy    = 0.0;
    TotalThermEnergy  = 0.0;
    TotalKinEnergy    = 0.0;
    TotalAllEnergy    = 0.0;

    [ Time, X1, X2, X3, CD, CS1, CS2, CS3, CE, CNe, CB1, CB2, CB3, CChi] = ...
        ReadMagnetofluidFields_Conserved( ProblemName, k, Directory);
    [ ~, ~, ~, ~, PD, PV1, PV2, PV3, PE, PNe, PB1, PB2, PB3, PChi] = ...
        ReadMagnetofluidFields_Primitive( ProblemName, k, Directory);
    [ ~, ~, ~, ~, AP, AT, AYe, AS, AE, AMe, AMp, AMn, AYp, AYn, AYa, AYh, AGm, ACs ] = ...
        ReadMagnetofluidFields_Auxiliary( ProblemName, k, Directory);
    [ ~, ~, ~, ~, DTCI, DShockX1, DShockX2, DShockX3, DTheta1, DTheta2, DTheta3, DMinE, DMaxE, DDiv ]...
        = ReadMagnetofluidFields_Diagnostic( ProblemName, k, Directory );   

    W = 1.0 ./ sqrt(1.0 - (PV1.^2 + PV2.^2 + PV3.^2));
    bu0 = W .* (CB1 .* PV1 + CB2 .* PV2 + CB3 .* PV3);
    bd0 = -bu0;
    PBSqr = bu0 .* bd0 + PB1.^2 + PB2.^2 + PB3.^2;
    
    MagEnergy = (W.^2 - 1.0/2.0) .* PBSqr - (bu0).^2;
    Div    = DDiv;
    DivEnergy   = 0.5 * CChi.^2;
    ThermEnergy = PE;
    KinEnergy   = PD .* W.^2 + AGm .* PE .* W.^2;
    AllEnergy   = MagEnergy + DivEnergy + ThermEnergy;
    
    MagEnergyAvg     = GetCellAverages(Dim, nNodes, MagEnergy);
    DivAvg           = abs(GetCellAverages(Dim, nNodes, Div));
    DivEnergyAvg     = GetCellAverages(Dim, nNodes, DivEnergy);
    ThermEnergyAvg   = GetCellAverages(Dim, nNodes, ThermEnergy);
    KinEnergyAvg     = GetCellAverages(Dim, nNodes, KinEnergy);
    AllEnergyAvg     = GetCellAverages(Dim, nNodes, AllEnergy);
    
    TotalMagEnergy   = sum(MagEnergyAvg   * dV, "all");
    TotalDiv         = sum(DivAvg         * dV, "all");
    TotalDivEnergy   = sum(DivEnergyAvg   * dV, "all");
    TotalThermEnergy = sum(ThermEnergyAvg * dV, "all");
    TotalKinEnergy   = sum(KinEnergyAvg   * dV, "all");
    TotalAllEnergy   = sum(AllEnergyAvg   * dV, "all");
    
    TotalMagEnergyVec(k+1)   = TotalMagEnergy;
    TotalDivVec(k+1)         = TotalDiv;
    TotalDivEnergyVec(k+1)   = TotalDivEnergy;
    TotalThermEnergyVec(k+1) = TotalThermEnergy;
    TotalKinEnergyVec(k+1)   = TotalKinEnergy;
    TotalAllEnergyVec(k+1)   = TotalAllEnergy;
    
    TimeVec(k+1) = Time;

end
figure(4)
semilogy(TimeVec, TotalMagEnergyVec,   'linewidth', 2.5);
hold on
semilogy(TimeVec, TotalDivEnergyVec,   'linewidth', 2.5);
semilogy(TimeVec, TotalThermEnergyVec, 'linewidth', 2.5);

figure(2)
plot(TimeVec, TotalMagEnergyVec - TotalMagEnergyVec(1), 'linewidth', 2.5)
hold on
plot(TimeVec, TotalThermEnergyVec - TotalThermEnergyVec(1), 'linewidth', 2.5)
plot(TimeVec, TotalDivEnergyVec - TotalDivEnergyVec(1), 'linewidth', 2.5)
plot(TimeVec, TotalKinEnergyVec - TotalKinEnergyVec(1), 'linewidth', 2.5)
plot(TimeVec, TotalKinEnergyVec + TotalMagEnergyVec - TotalKinEnergyVec(1) - TotalMagEnergyVec(1), 'linewidth', 2.5)

figure(3)
plot(TimeVec, TotalDivVec)

figure(5)
plot(TimeVec, TotalDivEnergyVec)