function [ E, T, Eta, J_I0, J_II0, J_I1, J_II1, OS] = readPairOpacityTable( opacityTableName )

  % Reads HDF5 Opacity Table

  disp( fprintf( 'INFO: Reading Pair Opacity from file: %s', opacityTableName ) );

  % Independent Variables:
  E   = h5read( opacityTableName, '/EnergyGrid/Values' );
  T   = h5read( opacityTableName, '/ThermoState/Temperature' );
  Eta = h5read( opacityTableName, '/EtaGrid/Values' );
  
  OS = h5read( opacityTableName, '/Scat_Pair_Kernels/Offsets' );
  
  Tmp = h5read( opacityTableName, '/Scat_Pair_Kernels/Kernels' );
  J_I0 = 10.^( squeeze(Tmp(:,:,1,:,:)) ) - OS(1);
  J_II0 = 10.^( squeeze(Tmp(:,:,2,:,:)) ) - OS(2);
  J_I1 = 10.^( squeeze(Tmp(:,:,3,:,:)) ) - OS(3);
  J_II1 = 10.^( squeeze(Tmp(:,:,4,:,:)) ) - OS(4);
  
end

