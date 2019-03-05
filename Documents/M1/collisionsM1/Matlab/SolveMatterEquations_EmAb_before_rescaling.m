function SolveMatterEquations_EmAb_before_rescaling( J, Chi, D, T, Y, E )

    % E_N       : nE_G x 1 (energy grid)
    % J, J_0    : nE_G x 1 (size of energy grid)
    % Chi       : nE_G x 1 (size of energy grid)
    % W2_N, W3_N: nE_G x 1 (size of energy grid)
    
    % constants needed
    global AtomicMassUnit BoltzmannConstant;
    
    iY = 1;
    iE = 2;
    U = zeros(2,1);
    C = zeros(2,1);
    FVEC = zeros(2,1);
    
    % \rho / m_B
    N_B = D / AtomicMassUnit;
    
    
    % (scaled) weights
    Theta2_N = 4 * pi * W2_N .* Chi ./ ( 1 + Chi );
    Theta3_N = 4 * pi * W3_N .* Chi ./ ( 1 + Chi );

    
    % compute chemical potential and derivatives
    [Mnu, dMnudT, dMnudY] = ComputeNeutrinoChemicalPotentials(D, T, Y);
    
    % equilibrium distribution
    J0 = FermiDirac( E_N, Mnu, BoltzmannConstant * T );

    % initial guess
    U(iY) = Y; U(iE) = E;
    
    % old states
    C(iY) = Theta2_N' * J + N_B * U(iY);
    C(iE) = Theta3_N' * J + N_B * U(iE);

   
    % --- Electron Fraction Equation ---

    FVEC(iY) = Theta2_N' * J0 + N_B * U(iY) - C(iY);

    % --- Internal Energy Equation ---

    FVEC(iE) = Theta3_N' * J0 + N_B * U(iE) - C(iE);

    % --- Scale Equations and Save Initial Evaluation ---

    FVEC = FVEC ./ C; 
    FVEC0 = FVEC;
    
    
    k = 0;
    CONVERGED = FALSE;
    while(~CONVERGED)

      k = k + 1;

      % QUESTION: The first output of ComputeSpecificInternalEnergy_TABLE 
      % should be the same as the input E? (i.e., is D T Y E in the input
      % self-consistent?)
      
      % compute derivatives
      [~, dEdT, dEdY] = ComputeSpecificInternalEnergy_TABLE(D, T, Y);
      
      % --- Derivative of J0 wrt. T (Constant Y) ---

      dJ0dT_Y = dFermiDiracdT( E_N, Mnu, BoltzmannConstant * T, dMnudT, T );

      % --- Derivative of J0 wrt. T (Constant T) ---

      dJ0dY_T = dFermiDiracdY( E_N, Mnu, BoltzmannConstant * T, dMnudY);

      % --- Derivative of J0 wrt. E (Constant Y) ---

      dJ0dE_Y = dJ0dT_Y / dEdT;

      % --- Derivative of J0 wrt. Y (Constant E) ---

      dJ0dY_E = dJ0dY_T - dJ0dT_Y * dEdY / dEdT;

      % --- Jacobian ---

      FJAC(1,1) = Theta2_N' * dJ0dY_E  + N_B;

      FJAC(1,2) = Theta2_N' * dJ0dE_Y;

      FJAC(2,1) = Theta3_N' * dJ0dY_E;

      FJAC(2,2) = Theta3_N' * dJ0dE_Y + N_B;

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

      Y = U(iY); E = U(iE);

      
      T = ComputeTemperatureFromSpecificInternalEnergy_TABLE(D, E, Y); 
            
      % --- Neutrino Chemical Potential and Derivatives ---

      [Mnu, dMnudT, dMnudY] = ComputeNeutrinoChemicalPotentials(D, T, Y);

      % --- Equilibrium Distribution ---

      J0 = FermiDirac( E_N, Mnu, BoltzmannConstant * T );

      % --- Electron Fraction Equation ---

      FVEC(iY) = Theta2_N' * J0 + N_B * U(iY) - C(iY);

      % --- Internal Energy Equation ---

      FVEC(iE) = Theta3_N' * J0 + N_B * U(iE) - C(iE);

      % --- Scale Equations ---

      FVEC = FVEC ./ C; 

      % --- Check for Convergence ---

      if((norm(FVEC) < Rtol * norm(FVEC0)) || (norm(dU./U) < Utol))

        CONVERGED = true;

        Iterations_Min = MIN( Iterations_Min, k );
        Iterations_Max = MAX( Iterations_Max, k );
        Iterations_Ave = Iterations_Ave + k;

      end

    

    end


end












