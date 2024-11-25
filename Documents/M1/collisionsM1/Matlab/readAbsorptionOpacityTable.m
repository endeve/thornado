function [ E, D, T, Y, Chi_N, Chi_AN, OS ] = readAbsorptionOpacityTable( opacityTableName )

  % Reads HDF5 Opacity Table
  % OS(1): offset for Neutrino
  % OS(2): offset for Antineutrino

  disp( fprintf( 'INFO: Reading Absorption Opacity from file: %s', opacityTableName ) );

  % Independent Variables:
  E = h5read( opacityTableName, '/EnergyGrid/Values' );
  D = h5read( opacityTableName, '/ThermoState/Density' );
  T = h5read( opacityTableName, '/ThermoState/Temperature' );
  Y = h5read( opacityTableName, '/ThermoState/Electron Fraction' );
  
  OS = h5read( opacityTableName, '/EmAb_CorrectedAbsorption/Offsets' );
  
  Tmp = h5read( opacityTableName, '/EmAb_CorrectedAbsorption/Electron Neutrino' );
  Chi_N = 10.^(Tmp) - OS(1); 
  
  Tmp = h5read( opacityTableName, '/EmAb_CorrectedAbsorption/Electron Antineutrino' );
  Chi_AN = 10.^(Tmp) - OS(2); 
  
end

