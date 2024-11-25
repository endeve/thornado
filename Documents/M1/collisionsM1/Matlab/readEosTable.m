function [ D, T, Y, nD, nT, nY, P, S, E, Me, Mp, Mn, Xp, Xn, Xa, Xh,...
           Zh, Ah, Eh, Eth, Gm, OS ] = readEosTable( eosTableName )

    % Reads HDF5 EOS Table

    disp( fprintf( 'INFO: Reading EOS from file: %s', eosTableName ) );

    % Independent Variables:
    D = h5read( eosTableName, '/ThermoState/Density' );
    T = h5read( eosTableName, '/ThermoState/Temperature' );
    Y = h5read( eosTableName, '/ThermoState/Electron Fraction' );
    
    Dims = h5read( eosTableName, '/ThermoState/Dimensions' );
    
    nD = Dims( 1 );
    nT = Dims( 2 );
    nY = Dims( 3 );

    minD = min( D );
    maxD = max( D );
    minT = min( T );
    maxT = max( T );
    minY = min( Y );
    maxY = max( Y );
    
    disp( fprintf( '  INFO: nD, nT, nY = %i, %i, %i', nD, nT, nY ) );
    
    disp( fprintf( '  INFO: Range in D: %d to %d', minD, maxD ) );
    disp( fprintf( '  INFO: Range in T: %d to %d', minT, maxT ) );
    disp( fprintf( '  INFO: Range in Y: %d to %d', minY, maxY ) );

    % Dependent Variables:
    
    OS = h5read( eosTableName, '/DependentVariables/Offsets' );
    
    % Pressure:
    iP = h5read( eosTableName, '/DependentVariables/iPressure' );
    P  = h5read( eosTableName, '/DependentVariables/Pressure' );
    P  = 10.^( P ) - OS(iP);
    
    % Entropy Per Baryon:
    iS = h5read( eosTableName, '/DependentVariables/iEntropyPerBaryon' );
    S  = h5read( eosTableName, '/DependentVariables/Entropy Per Baryon' );
    S  = 10.^( S ) - OS(iS);
    
    % Ineternal Energy Density:
    E.i = h5read( eosTableName, '/DependentVariables/iInternalEnergyDensity' );
    E.v3D  = h5read( eosTableName, '/DependentVariables/Internal Energy Density' );
    E.v3D  = 10.^( E.v3D ) - OS(E.i);
    
    % Electron Chemical Potential:
    Me.i = h5read( eosTableName, '/DependentVariables/iElectronChemicalPotential' );
    Me.v3D  = h5read( eosTableName, '/DependentVariables/Electron Chemical Potential' );
    Me.v3D  = 10.^( Me.v3D ) - OS(Me.i);
    
    % Proton Chemical Potential:
    Mp.i = h5read( eosTableName, '/DependentVariables/iProtonChemicalPotential' );
    Mp.v3D  = h5read( eosTableName, '/DependentVariables/Proton Chemical Potential' );
    Mp.v3D  = 10.^( Mp.v3D ) - OS(Mp.i);
    
    % Neutron Chemical Potential:
    Mn.i = h5read( eosTableName, '/DependentVariables/iNeutronChemicalPotential' );
    Mn.v3D  = h5read( eosTableName, '/DependentVariables/Neutron Chemical Potential' );
    Mn.v3D  = 10.^( Mn.v3D ) - OS(Mn.i);
    
    % Proton Mass Fraction:
    iXp = h5read( eosTableName, '/DependentVariables/iProtonMassFraction' );
    Xp  = h5read( eosTableName, '/DependentVariables/Proton Mass Fraction' );
    Xp  = 10.^( Xp ) - OS(iXp);
    
    % Neutron Mass Fraction:
    iXn = h5read( eosTableName, '/DependentVariables/iNeutronMassFraction' );
    Xn  = h5read( eosTableName, '/DependentVariables/Neutron Mass Fraction' );
    Xn  = 10.^( Xn ) - OS(iXn);
    
    % Alpha Mass Fraction:
    iXa = h5read( eosTableName, '/DependentVariables/iAlphaMassFraction' );
    Xa  = h5read( eosTableName, '/DependentVariables/Alpha Mass Fraction' );
    Xa  = 10.^( Xa ) - OS(iXa);
    
    % Heavy Mass Fraction:
    iXh = h5read( eosTableName, '/DependentVariables/iHeavyMassFraction' );
    Xh  = h5read( eosTableName, '/DependentVariables/Heavy Mass Fraction' );
    Xh  = 10.^( Xh ) - OS(iXh);
    
    % Heavy Charge Number:
    iZh = h5read( eosTableName, '/DependentVariables/iHeavyChargeNumber' );
    Zh  = h5read( eosTableName, '/DependentVariables/Heavy Charge Number' );
    Zh  = 10.^( Zh ) - OS(iZh);
    
    % Heavy Mass Number:
    iAh = h5read( eosTableName, '/DependentVariables/iHeavyMassNumber' );
    Ah  = h5read( eosTableName, '/DependentVariables/Heavy Mass Number' );
    Ah  = 10.^( Ah ) - OS(iAh);
    
    % Heavy Binding Energy:
    iEh = h5read( eosTableName, '/DependentVariables/iHeavyBindingEnergy' );
    Eh  = h5read( eosTableName, '/DependentVariables/Heavy Binding Energy' );
    Eh  = 10.^( Eh ) - OS(iEh);
    
    % Thermal Energy:
    iEth = h5read( eosTableName, '/DependentVariables/iThermalEnergy' );
    Eth  = h5read( eosTableName, '/DependentVariables/Thermal Energy' );
    Eth  = 10.^( Eth ) - OS(iEth);
    
    % Gamma1:
    iGm = h5read( eosTableName, '/DependentVariables/iGamma1' );
    Gm  = h5read( eosTableName, '/DependentVariables/Gamma1' );
    Gm  = 10.^( Gm ) - OS(iGm);
    
end

