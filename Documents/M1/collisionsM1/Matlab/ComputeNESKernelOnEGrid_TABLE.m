function  [Chi] = ComputeNESKernelOnEGrid_TABLE(E_N, D, T, Y) 
    
    global g_D1D g_T1D g_Y1D g_E1D g_Chi_N g_Ab_OS;
    
    [ Me, ~, ~, ~ ] = interpolateDifferentiateEos...
        ( D, T, Y, g_D1D, g_T1D, g_Y1D, g_Me.v3D, g_EOS_OS(g_Me.i) );
    
    [ Chi ] = interpolate4D(E_N, D, T, Y, g_E1D, g_D1D, ...
        g_T1D, g_Y1D, [1 1 1 0], g_Chi_N, g_Ab_OS(1));

end