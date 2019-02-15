function readTables()
% read tables
    
    % grids
    global g_D1D g_T1D g_Y1D;
    
    % chemical potential
    global g_Me g_Mp g_Mn;
    
    % internal energy density
    global g_E;
    
    % EOS offset
    global g_EOS_OS;
    
    % Absorption (neutrino and antineutrino)
    global g_Chi_N g_Chi_AN g_Ab_OS;
    
    % energy grid
    global g_E1D;
    
    [ g_D1D, g_T1D, g_Y1D, nD, nT, nY, P, S, g_E, g_Me, g_Mp, g_Mn, Xp, Xn, Xa, Xh,...
    Zh, Ah, Eh, Eth, Gm, g_EOS_OS ] = readEosTable('wl-EOS-SFHo-15-25-50.h5');

    [ g_E1D, ~, ~, ~, g_Chi_N, g_Chi_AN, g_Ab_OS ] = readAbsorptionOpacityTable('wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5');
    
     
% [ E, D, T, Y, R_0_N, R_1_N, R_0_AN, R_1_AN, OS ] = readIsoScatteringOpacityTable('wl-Op-SFHo-15-25-50-E40-B85-Iso.h5');

% [ E, T, Eta, H_I0, H_II0, H_I1, H_II1, OS] = readNesScatteringOpacityTable('wl-Op-SFHo-15-25-50-E40-B85-NES.h5');

end