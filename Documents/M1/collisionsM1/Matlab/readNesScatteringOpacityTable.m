function [ E, T, Eta, H_I0, H_II0, H_I1, H_II1, OS] = readNesScatteringOpacityTable( opacityTableName )

  % Reads HDF5 Opacity Table

  disp( fprintf( 'INFO: Reading NES Scattering Opacity from file: %s', opacityTableName ) );

  % Independent Variables:
  E   = h5read( opacityTableName, '/EnergyGrid/Values' );
  T   = h5read( opacityTableName, '/ThermoState/Temperature' );
  Eta = h5read( opacityTableName, '/EtaGrid/Values' );
  
  OS = h5read( opacityTableName, '/Scat_NES_Kernels/Offsets' );
  
  Tmp = h5read( opacityTableName, '/Scat_NES_Kernels/Kernels' );
  H_I0 = 10.^( squeeze(Tmp(:,:,1,:,:)) ) - OS(1);
  H_II0 = 10.^( squeeze(Tmp(:,:,2,:,:)) ) - OS(2);
  H_I1 = 10.^( squeeze(Tmp(:,:,3,:,:)) ) - OS(3);
  H_II1 = 10.^( squeeze(Tmp(:,:,4,:,:)) ) - OS(4);
  
end

