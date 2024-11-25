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
    
    % Emission and Absorption (neutrino and antineutrino)
    global g_Chi_Ne g_Chi_ANe g_Ab_OS;
    
    % energy grid
    global g_E1D;
    
    % Eta grid
    global g_Eta;
    
    % NES kernels
    global g_H_I0 g_H_II0 g_H_I1 g_H_II1; 
    
    % NES offset
    global g_NES_OS;
    
    % NES kernels
    global g_J_I0 g_J_II0 g_J_I1 g_J_II1; 
    
    % NES offset
    global g_Pair_OS;

    [ g_D1D, g_T1D, g_Y1D, nD, nT, nY, P, S, g_E, g_Me, g_Mp, g_Mn, Xp, Xn, Xa, Xh,...
    Zh, Ah, Eh, Eth, Gm, g_EOS_OS ] = readEosTable('wl-EOS-SFHo-15-25-50.h5');

    [ g_E1D, ~, ~, ~, g_Chi_Ne, g_Chi_ANe, g_Ab_OS ] = readAbsorptionOpacityTable('wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5');
    
     
% [ E, D, T, Y, R_0_N, R_1_N, R_0_AN, R_1_AN, OS ] = readIsoScatteringOpacityTable('wl-Op-SFHo-15-25-50-E40-B85-Iso.h5');

    [ ~, ~, g_Eta, g_H_I0, g_H_II0, g_H_I1, g_H_II1, g_NES_OS] = readNesScatteringOpacityTable('wl-Op-SFHo-15-25-50-E40-B85-NES.h5');

    [ ~, ~, ~, g_J_I0, g_J_II0, g_J_I1, g_J_II1, g_Pair_OS] = readPairOpacityTable('wl-Op-SFHo-15-25-50-E40-B85-Pair.h5');

end