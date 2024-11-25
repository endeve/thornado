function  [R_Pr_Ne, R_An_Ne, R_Pr_ANe, R_An_ANe] = ComputePairOpacityOnEGrid_TABLE(E_N, D, T, Y) 
    
    global g_D1D g_T1D g_Y1D g_E1D g_Me g_EOS_OS g_Eta g_J_I0 g_J_II0 g_Pair_OS;
%     global g_H_I1 g_H_II1;
    global BoltzmannConstant;
    
    kT = BoltzmannConstant * T;
    
    [ Me, ~, ~, ~ ] = interpolateDifferentiateEos...
        ( D, T, Y, g_D1D, g_T1D, g_Y1D, g_Me.v3D, g_EOS_OS(g_Me.i) );
    
    Eta = Me / kT;
    
    % zero-th order for now
    [ J_I0 ] = interpolate2D2D( E_N, E_N, T, Eta, g_E1D, g_E1D, g_T1D, g_Eta, [1 1 1 1], g_J_I0, g_Pair_OS(1));
    [ J_II0 ] = interpolate2D2D( E_N, E_N, T, Eta, g_E1D, g_E1D, g_T1D, g_Eta, [1 1 1 1], g_J_II0, g_Pair_OS(2));
%     [ J_I1 ] = interpolate2D2D( E_N, E_N, T, Eta, g_E1D, g_E1D, g_T1D, g_Eta, [1 1 1 1], g_J_I1, g_Pair_OS(3));
%     [ J_II1 ] = interpolate2D2D( E_N, E_N, T, Eta, g_E1D, g_E1D, g_T1D, g_Eta, [1 1 1 1], g_J_II1, g_Pair_OS(4));


    % Constants from (C50) in Bruenn 1985.
    % see also wlOpacityInterpolationModule.f90
    Cv = 0.96; Ca = 0.5;
    
    % this is (77) in the Tech Report, check the scaling
    R_An_Ne = (Cv + Ca)^2 * J_I0 + (Cv - Ca)^2 * J_II0;
    R_An_ANe = (Cv - Ca)^2 * J_I0 + (Cv + Ca)^2 * J_II0;
 
%%%
% Potential saving: interpolating only upper half, computing the upper half 
% of R_An_Ne and  R_An_ANe, and getting the lower half of R_An_Ne and 
% R_An_ANe from R_An_ANe' and R_An_Ne', respectively.
% Can do the same for R_Pr
%%%
    
    for iE  = 1:length(E_N)
        for iEp = iE+1:length(E_N)
            % assign the lower triangular part to enforce symmetry
            R_An_Ne(iEp,iE) = (Cv + Ca)^2 * J_II0(iE,iEp) + (Cv - Ca)^2 * J_I0(iE,iEp);
            
            R_An_ANe(iEp,iE) = (Cv - Ca)^2 * J_II0(iE,iEp) + (Cv + Ca)^2 * J_I0(iE,iEp);
        end
    end
    
    R_Pr_Ne = zeros(size(R_An_Ne));
    R_Pr_ANe = zeros(size(R_An_ANe));
    

%%%
% Potential saving: Precompute the exp part and save as a matrix.
%%%
    for iE  = 1:length(E_N)
        for iEp = 1:length(E_N)
            R_Pr_Ne(iE,iEp) = R_An_Ne(iE,iEp) * exp(-(E_N(iE) + E_N(iEp)) / kT);
            
            R_Pr_ANe(iE,iEp) = R_An_ANe(iE,iEp) * exp(-(E_N(iE) + E_N(iEp)) / kT);
        end
    end
%     R_Pr_ANe = R_Pr_Ne';
%     R_An_ANe = R_An_Ne';

end