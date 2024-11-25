function [M, dMdT, dMdY] = ComputeNeutrinoChemicalPotentials(D, T, Y)

    global g_D1D g_T1D g_Y1D g_Me g_Mp g_Mn g_EOS_OS;

    % use D T Y to find Me Mp Mn
    [ Me, ~, dMedT, dMedY ] = interpolateDifferentiateEos...
        ( D, T, Y, g_D1D, g_T1D, g_Y1D, g_Me.v3D, g_EOS_OS(g_Me.i) );
    [ Mp, ~, dMpdT, dMpdY ] = interpolateDifferentiateEos...
        ( D, T, Y, g_D1D, g_T1D, g_Y1D, g_Mp.v3D, g_EOS_OS(g_Mp.i) );
    [ Mn, ~, dMndT, dMndY ] = interpolateDifferentiateEos...
        ( D, T, Y, g_D1D, g_T1D, g_Y1D, g_Mn.v3D, g_EOS_OS(g_Mn.i) );

    M = ( Me + Mp ) - Mn;
    dMdT = ( dMedT + dMpdT ) - dMndT;
    dMdY = ( dMedY + dMpdY ) - dMndY;

end