function [ E, D, T, Y, R_0_N, R_1_N, R_0_AN, R_1_AN, OS ] = readIsoScatteringOpacityTable( opacityTableName )

  % Reads HDF5 Opacity Table

  disp( fprintf( 'INFO: Reading Isoenergetic Opacity from file: %s', opacityTableName ) );

  % Independent Variables:
  E = h5read( opacityTableName, '/EnergyGrid/Values' );
  D = h5read( opacityTableName, '/ThermoState/Density' );
  T = h5read( opacityTableName, '/ThermoState/Temperature' );
  Y = h5read( opacityTableName, '/ThermoState/Electron Fraction' );
  
  OS = h5read( opacityTableName, '/Scat_Iso_Kernels/Offsets' );
  OS = reshape(OS, [4,1]);
  Tmp = h5read( opacityTableName, '/Scat_Iso_Kernels/Electron Neutrino' );
  R_0_N = 10.^( squeeze(Tmp(:,1,:,:,:)) ) -OS(1);
  R_1_N = 10.^( squeeze(Tmp(:,2,:,:,:)) ) -OS(3);

  Tmp = h5read( opacityTableName, '/Scat_Iso_Kernels/Electron Antineutrino' );
  R_0_AN = 10.^( squeeze(Tmp(:,1,:,:,:)) ) -OS(2);
  R_1_AN = 10.^( squeeze(Tmp(:,2,:,:,:)) ) -OS(4);
  
end

