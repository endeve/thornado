function readTestConstants()

global AtomicMassUnit PlanckConstant BoltzmannConstant SpeedOfLight;

global Erg2MeV

% AtomicMassUnit = 1.0;
% PlanckConstant = 1.0;
% BoltzmannConstant = 1.0;
% SpeedOfLight = 1.0;


SpeedOfLightCGS          = 2.99792458e+10;   % cm / sec
BoltzmannConstantCGS     = 1.3806503e-16;   % erg / K
PlanckConstantCGS        = 6.62606876e-27;  % erg * s

% MKS
ElectronVoltMKS          = 1.602176462e-19;
% BoltzmannConstantMKS     = 1.3806503e-23;
% PlanckConstantMKS        = 6.62606876e-34;
AvogadroConstantMKS      = 6.02214199e23;


Erg2MeV = 1/(ElectronVoltMKS * 1e+6 * 1e+7);

SpeedOfLight        = SpeedOfLightCGS;
PlanckConstant      = PlanckConstantCGS * Erg2MeV;           % MeV * K
BoltzmannConstant   = BoltzmannConstantCGS * Erg2MeV;     % MeV / K
% BoltzmannConstant   = 8.61733261e-11;   
AtomicMassUnit      = 1.0 / AvogadroConstantMKS;

% GravitationalConstant = 1.0;

 
% --- Length ---
% Centimeter = 1.0;
% Meter      = 1.0e+2 * Centimeter;
% Kilometer  = 1.0e+3 * Meter;
% 

% --- Time ---
% Second      = SpeedOfLightCGS / SpeedOfLight * Centimeter;
% Millisecond = 1.0e-3 * Second;
% Microsecond = 1.0e-6 * Second;
% 

% --- Mass ---
% Gram        = 1.0; % Double-check this, do we have to use GravitationalConstant?
% Kilogram    = 1.0e+3 * Gram;
% 

% Kilogram  = GravitationalConstantMKS / GravitationalConstant * Meter^3 / Second^2;
% Gram      = 1.0e-3 * Kilogram;
% SolarMass = 1.98892e30 * Kilogram;
% 


% % --- Other Units of Measure and Constants ---
% Erg            = Gram * ( Centimeter / Second )^2;
% Joule          = Kilogram * ( Meter / Second )^2;
% ElectronVolt   = ElectronVoltMKS * Joule;
% MeV            = 1.0e6 * ElectronVolt;
% Kelvin         = BoltzmannConstantMKS / BoltzmannConstant * Joule;
% PlanckConstant = PlanckConstantMKS * Joule * Second;
% AtomicMassUnit = Gram / AvogadroConstantMKS;

% 
% Joule          = Kilogram * ( Meter / Second )^2;
% Erg            = Gram * ( Centimeter / Second )^2;
% Bethe          = 1.0e51 * Erg;
% ElectronVolt   = ElectronVoltMKS * Joule;
% MeV            = 1.0e6 * ElectronVolt;
% Kelvin         = BoltzmannConstantMKS / BoltzmannConstant * Joule;
% Newton         = Joule / Meter;
% Dyne           = Erg / Centimeter;
% PlanckConstant = PlanckConstantMKS * Joule * Second;
% AtomicMassUnit = Gram / AvogadroConstantMKS;


end


% --- MKS ---
% SpeedOfLightMKS          = 2.99792458e8;
% GravitationalConstantMKS = 6.673e-11;
% BoltzmannConstantMKS     = 1.3806503e-23;
% ElectronVoltMKS          = 1.602176462e-19;
% PlanckConstantMKS        = 6.62606876e-34;
% AvogadroConstantMKS      = 6.02214199e23;
% 