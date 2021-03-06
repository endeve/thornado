function [J0, D, T, Y, E, iter] = SolveMatterEquations_EmAb( J, Chi, D, T, Y, E )

    % g_E_N       : nE_G x 1 (energy grid)
    % J, J_0    : nE_G x 1 (size of energy grid)
    % Chi       : nE_G x 1 (size of energy grid)
    % W2_N, W3_N: nE_G x 1 (size of energy grid)
    
    % constants
    global AtomicMassUnit BoltzmannConstant Erg2MeV SpeedOfLight PlanckConstant;

    % energy grid
    global g_E_N g_W2_N g_W3_N;
    
    % iteration data
    
    global g_Iterations_Min g_Iterations_Max g_Iterations_Ave;
    
    hc = PlanckConstant * SpeedOfLight;
    c = SpeedOfLight;
    
    Rtol = 1e-8; Utol = 1e-10; maxIter = 100;

    iY = 1;
    iE = 2;
    U = zeros(2,1);
    C = zeros(2,1);
    FVEC = zeros(2,1);
    
       
    % \rho / m_B
    N_B = D / AtomicMassUnit;
    
    % scales
    s_Y = N_B;
    s_E = D;
    
    Chi = Chi * c;

    
%     % weights
%     Numbins = 20;
%     NumPts = 2;
%     [E_N, W2_N, W3_N] = ComputePointsAndWeightsE(Numbins,NumPts);
    
    % (scaled) weights
    Theta2_N = 4 * pi * g_W2_N .* Chi ./ ( 1 + Chi );
    Theta3_N = 4 * pi * g_W3_N .* Chi ./ ( 1 + Chi );

    Theta2_N = Theta2_N / (hc)^3;
    Theta3_N = Theta3_N / (hc)^3;



    
    % compute chemical potential and derivatives
    [Mnu, dMnudT, dMnudY] = ComputeNeutrinoChemicalPotentials(D, T, Y);
    
    % equilibrium distribution
    J0 = FermiDirac( g_E_N, Mnu, BoltzmannConstant * T );

    % store initial
    Y0 = Y; 
%     E0 = E;
    E0 = E * Erg2MeV; % change unit

    
    % scaling for Jacobian
    s_Y = s_Y * Y0;
    s_E = s_E * E0;

    
    % initial guess
    U(iY) = Y / Y0; 
    U(iE) = E * Erg2MeV / E0; % change unit
    
    % old states
    C(iY) = Theta2_N' * J + s_Y * U(iY);
    C(iE) = Theta3_N' * J + s_E * U(iE);

   
    % --- Electron Fraction Equation ---

    FVEC(iY) = Theta2_N' * J0 + s_Y * U(iY) - C(iY);

    % --- Internal Energy Equation ---

    FVEC(iE) = Theta3_N' * J0 + s_E * U(iE) - C(iE);

    % --- Scale Equations and Save Initial Evaluation ---

    FVEC = FVEC ./ C; 
    FVEC0 = FVEC;
    
    
    k = 0;
    CONVERGED = false;
    while(~CONVERGED)

      k = k + 1;

      % QUESTION: The first output of ComputeSpecificInternalEnergy_TABLE 
      % should be the same as the input E? (i.e., is D T Y E in the input
      % self-consistent?) --> YES!
      
      % compute derivatives
      [~, dEdT, dEdY] = ComputeSpecificInternalEnergy_TABLE(D, T, Y);
      
      % change units
      dEdT = dEdT * Erg2MeV / E0 ;
      dEdY = dEdY * Y0 * Erg2MeV / E0;
      
      % --- Derivative of J0 wrt. T (Constant Y) ---

      dJ0dT_Y = dFermiDiracdT(g_E_N, Mnu, BoltzmannConstant * T, dMnudT, T);

      % --- Derivative of J0 wrt. Y (Constant T) ---

      dJ0dY_T = dFermiDiracdY(g_E_N, Mnu, BoltzmannConstant * T, dMnudY);

      dJ0dY_T = dJ0dY_T * Y0;
      
      % --- Derivative of J0 wrt. E (Constant Y) ---

      dJ0dE_Y = dJ0dT_Y / dEdT;

      % --- Derivative of J0 wrt. Y (Constant E) ---

      dJ0dY_E = dJ0dY_T - dJ0dT_Y * dEdY / dEdT;

      % --- Jacobian ---

      FJAC(1,1) = Theta2_N' * dJ0dY_E  + s_Y;

      FJAC(1,2) = Theta2_N' * dJ0dE_Y;

      FJAC(2,1) = Theta3_N' * dJ0dY_E;

      FJAC(2,2) = Theta3_N' * dJ0dE_Y + s_E;

      % --- Scale Jacobian ---

      FJAC(:,1) = FJAC(:,1) ./ C;

      FJAC(:,2) = FJAC(:,2) ./ C;

      % --- Determinant of Jacobian ---

      DJAC = FJAC(1,1) * FJAC(2,2) - FJAC(2,1) * FJAC(1,2);

      % --- Invert Jacobian ---

      IJAC(1,1) =   FJAC(2,2) / DJAC;
      IJAC(2,1) = - FJAC(2,1) / DJAC;
      IJAC(1,2) = - FJAC(1,2) / DJAC;
      IJAC(2,2) =   FJAC(1,1) / DJAC;

      % --- Correction ---

      dU = - IJAC * FVEC;

      % --- Apply Correction ---

      U = U + dU;

      Y = U(iY) * Y0; 
      E = U(iE) * E0 / Erg2MeV; % change unit

      
      T = ComputeTemperatureFromSpecificInternalEnergy_TABLE(D, E, Y); 
            
      % --- Neutrino Chemical Potential and Derivatives ---

      [Mnu, dMnudT, dMnudY] = ComputeNeutrinoChemicalPotentials(D, T, Y);

      % --- Equilibrium Distribution ---

      J0 = FermiDirac( g_E_N, Mnu, BoltzmannConstant * T );

      % --- Electron Fraction Equation ---

      FVEC(iY) = Theta2_N' * J0 + s_Y * U(iY) - C(iY);

      % --- Internal Energy Equation ---

      FVEC(iE) = Theta3_N' * J0 + s_E * U(iE) - C(iE);

      % --- Scale Equations ---

      FVEC = FVEC ./ C; 

      % --- Check for Convergence ---

      if((norm(FVEC) < Rtol * norm(FVEC0)) || (norm(dU./U) < Utol))
%       if((norm(FVEC) < Rtol * norm(FVEC0)))

        CONVERGED = true;

        g_Iterations_Min = min( g_Iterations_Min, k );
        g_Iterations_Max = max( g_Iterations_Max, k );
        g_Iterations_Ave = g_Iterations_Ave + k;

      end
      if(k >= maxIter)
         disp("Failed to converge within maxIter.");
         break;
      end
    

    end
    iter = k;

end












