function  [E, dEdT, dEdY] = ComputeSpecificInternalEnergy_TABLE(D, T, Y) 
    
    global g_D1D g_T1D g_Y1D g_E g_EOS_OS;

    [ E, ~, dEdT, dEdY ] = interpolateDifferentiateEos...
        ( D, T, Y, g_D1D, g_T1D, g_Y1D, g_E.v3D, g_EOS_OS(g_E.i) );

end