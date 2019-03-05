function readConstants()

%global SpeedOfLightMKS GravitationalConstantMKS BoltzmannConstantMKS ElectronVoltMKS PlanckConstantMKS AvogadroConstantMKS

global AtomicMassUnit PlanckConstant BoltzmannConstant SpeedOfLight;

SpeedOfLightMKS          = 2.99792458e8;
GravitationalConstantMKS = 6.673e-11;
BoltzmannConstantMKS     = 1.3806503e-23;
ElectronVoltMKS          = 1.602176462e-19;
PlanckConstantMKS        = 6.62606876e-34;
AvogadroConstantMKS      = 6.02214199e23;

SpeedOfLight          = 1.0;
GravitationalConstant = 1.0;
BoltzmannConstant     = 1.0;

% --- Length ---

Meter      = 1.0;
Centimeter = 1.0e-2 * Meter;
Kilometer  = 1.0e+3 * Meter;

% --- Time ---

Second      = SpeedOfLightMKS / SpeedOfLight * Meter;
Millisecond = 1.0e-3 * Second;
Microsecond = 1.0e-6 * Second;

% --- Mass ---

Kilogram  = GravitationalConstantMKS / GravitationalConstant * Meter^3 / Second^2;
Gram      = 1.0e-3 * Kilogram;
SolarMass = 1.98892e30 * Kilogram;

% --- Other Units of Measure and Constants ---

Joule          = Kilogram * ( Meter / Second )^2;
Erg            = Gram * ( Centimeter / Second )^2;
Bethe          = 1.0e51 * Erg;
ElectronVolt   = ElectronVoltMKS * Joule;
MeV            = 1.0e6 * ElectronVolt;
Kelvin         = BoltzmannConstantMKS / BoltzmannConstant * Joule;
Newton         = Joule / Meter;
Dyne           = Erg / Centimeter;
PlanckConstant = PlanckConstantMKS * Joule * Second;
AtomicMassUnit = Gram / AvogadroConstantMKS;

end
