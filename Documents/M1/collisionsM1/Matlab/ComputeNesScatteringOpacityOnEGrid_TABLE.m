function  [R_in_NES, R_out_NES] = ComputeNesScatteringOpacityOnEGrid_TABLE(E_N, D, T, Y) 
    
    global g_D1D g_T1D g_Y1D g_E1D g_Me g_EOS_OS g_Eta g_H_I0 g_H_II0 g_NES_OS;
%     global g_H_I1 g_H_II1;
    global BoltzmannConstant;
    
    kT = BoltzmannConstant * T;
    
    [ Me, ~, ~, ~ ] = interpolateDifferentiateEos...
        ( D, T, Y, g_D1D, g_T1D, g_Y1D, g_Me.v3D, g_EOS_OS(g_Me.i) );
    
    Eta = Me / kT;
    
    % zero-th order for now
    [ H_I0 ] = interpolate2D2D( E_N, E_N, T, Eta, g_E1D, g_E1D, g_T1D, g_Eta, [1 1 1 1], g_H_I0, g_NES_OS(1));
    [ H_II0 ] = interpolate2D2D( E_N, E_N, T, Eta, g_E1D, g_E1D, g_T1D, g_Eta, [1 1 1 1], g_H_II0, g_NES_OS(2));
%     [ H_I1 ] = interpolate2D2D( E_N, E_N, T, Eta, g_E1D, g_E1D, g_T1D, g_Eta, [1 1 1 1], g_H_I1, g_NES_OS(3));
%     [ H_II1 ] = interpolate2D2D( E_N, E_N, T, Eta, g_E1D, g_E1D, g_T1D, g_Eta, [1 1 1 1], g_H_II1, g_NES_OS(4));
    
    % not sure how to compute this. (C50) in Bruenn 1985.
    % check wlOpacityInterpolationModule.f90
    Cv = 0.96; Ca = 0.5;
    
    % this is (77) in the Tech Report, check the scaling
    R_out_NES = (Cv + Ca)^2 * H_I0 + (Cv - Ca)^2 * H_II0;
    
    
    for iE  = 1:length(E_N)
        for iEp = iE+1:length(E_N)
            % assign the lower triangular part to enforce detailed balance
            R_out_NES(iEp,iE) = R_out_NES(iE,iEp) * exp( ( E_N(iE) - E_N(iEp) ) / kT );
%             R_out_NES(iEp,iE) = R_out_NES(iE,iEp) * exp( -( E_N(iE) - E_N(iEp) ) / kT ); 
        end
    end
    R_in_NES = R_out_NES';

end