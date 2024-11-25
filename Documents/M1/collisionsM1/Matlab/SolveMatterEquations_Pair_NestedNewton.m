function [J, J0, D, T, Y, E, iter, iter_in] = SolveMatterEquations_Pair_NestedNewton( Jin, dt, Chi, D, T, Y, E)

% g_E_N       : nE_G x 1 (energy grid)
% J, J_0    : nE_G x 1 (size of energy grid)
% Chi       : nE_G x 1 (size of energy grid)
% W2_N, W3_N: nE_G x 1 (size of energy grid)

% constants
global AtomicMassUnit BoltzmannConstant Erg2MeV SpeedOfLight PlanckConstant;

% energy grid
global g_E_N g_W2_N g_W3_N;


hc = PlanckConstant * SpeedOfLight;
c = SpeedOfLight;

Rtol = 1e-8; Utol = 1e-10; maxIter = 100;
iter_in = 0;

iY = 1;
iE = 2;
iJNe = 1:length(Jin.Ne);
iJANe = iJNe(end) + (1:length(Jin.ANe));

U = zeros(2,1);
C = zeros(2,1);


% \rho / m_B
N_B = D / AtomicMassUnit;

Chi.Ne = Chi.Ne * c * dt;
Chi.ANe = Chi.ANe * c * dt;

% scales
s_Y = N_B;
s_E = D;

% weights (for NES)
W2_N = g_W2_N;


% (scaled) weights
Theta2_N = 4 * pi * g_W2_N;
Theta3_N = 4 * pi * g_W3_N;

Theta2_N = Theta2_N / (hc)^3;
Theta3_N = Theta3_N / (hc)^3;

% store initial values
Y0 = Y;
E0 = E * Erg2MeV; % change unit

% update scales
s_Y = s_Y * Y0;
s_E = s_E * E0;

% initial guess (both 1 for now)
U(iY) = Y / Y0;
U(iE) = E * Erg2MeV / E0; % change unit

U0 = U;

J = Jin;


% calculate known values
C(iY) = Theta2_N' * (Jin.Ne - Jin.ANe)/ s_Y;
C(iE) = Theta3_N' * (Jin.Ne + Jin.ANe) / s_E;


k = 0;
CONVERGED = false;

% Anderson acceleration truncation parameter
m = 3;
FVEC = zeros(length(J.Ne)+length(J.ANe),1);
FJAC = zeros(length(J.Ne)+length(J.ANe));
Uvec = zeros(length(U),m);
Fvec = zeros(length(U),m);
alpha = zeros(m,1);


% compute chemical potential and derivatives
[Mnu, ~, ~] = ComputeNeutrinoChemicalPotentials(D, T, Y);

% equilibrium distribution
[J0.Ne, J0.ANe] = FermiDirac( g_E_N, Mnu, BoltzmannConstant * T );
    
% FOR IMPLICIT R_NES, UPDATE (interpolate) R_NES here:
[R_in_NES.Ne, R_out_NES.Ne, R_in_NES.ANe, R_out_NES.ANe] = ComputeNesScatteringOpacityOnEGrid_TABLE(g_E_N, D, T, Y);

% scale R_NES by c * dt
R_in_NES.Ne = R_in_NES.Ne * c * dt;
R_out_NES.Ne = R_out_NES.Ne * c * dt;

R_in_NES.ANe = R_in_NES.ANe * c * dt;
R_out_NES.ANe = R_out_NES.ANe * c * dt;

eta_NES.Ne = R_in_NES.Ne' * (W2_N .* J.Ne);
Chi_NES.Ne = eta_NES.Ne + R_out_NES.Ne' * (W2_N .* (1 - J.Ne));

eta_NES.ANe = R_in_NES.ANe' * (W2_N .* J.ANe);
Chi_NES.ANe = eta_NES.ANe + R_out_NES.ANe' * (W2_N .* (1 - J.ANe));


% update Pair kernels (R_Pr and R_An)
[R_Pr.Ne, R_An.Ne, R_Pr.ANe, R_An.ANe] = ComputePairOpacityOnEGrid_TABLE(g_E_N, D, T, Y);

% scale R_Pr and R_An by c * dt
R_Pr.Ne = R_Pr.Ne * c * dt;
R_An.Ne = R_An.Ne * c * dt;

R_Pr.ANe = R_Pr.ANe * c * dt;
R_An.ANe = R_An.ANe * c * dt;

eta_TP.Ne = R_Pr.Ne' * (W2_N .* (1 - J.ANe));
Chi_TP.Ne = eta_TP.Ne + R_An.Ne' * (W2_N .* J.ANe);
% 
eta_TP.ANe = R_Pr.ANe' * (W2_N .* (1 - J.Ne));
Chi_TP.ANe = eta_TP.ANe + R_An.ANe' * (W2_N .* J.Ne);

% update J for electron type neutrino and antineutrino
% J.Ne = (Jin.Ne + Chi.Ne.*J0.Ne + eta_NES.Ne + eta_TP.Ne)./(1 + Chi.Ne + Chi_NES.Ne + Chi_TP.Ne);
% J.ANe = (Jin.ANe + Chi.ANe.*J0.ANe + eta_NES.ANe + eta_TP.ANe)./(1 + Chi.ANe + Chi_NES.ANe + Chi_TP.ANe);

    
while((~CONVERGED)&&(k<=maxIter))
    
    k = k + 1;
    mk = min(m,k);
 
    k_in = 0;
    CONVERGED_IN = false;
    J_in = J;
    
    % Solve inner loop
    while((~CONVERGED_IN)&&(k_in<=maxIter))
        
        % NEWTON
        
                
        k_in = k_in + 1;
        
        
        FVEC(iJNe) = J.Ne - Jin.Ne + Chi.Ne.*(J.Ne - J0.Ne) ...
            - ((1 - J.Ne).* (R_in_NES.Ne'*(W2_N .* J.Ne)) - J.Ne.*(R_out_NES.Ne'*(W2_N .* (1-J.Ne)))) ...
            - ((1 - J.Ne).* (R_Pr.Ne'*(W2_N .* (1 - J.ANe))) - J.Ne.*(R_An.Ne'*(W2_N .* J.ANe)));
        FVEC(iJANe) = J.ANe - Jin.ANe + Chi.ANe.*(J.ANe - J0.ANe) ...
            - ((1 - J.ANe).* (R_in_NES.ANe'*(W2_N .* J.ANe)) - J.ANe.*(R_out_NES.ANe'*(W2_N .* (1-J.ANe)))) ...
            - ((1 - J.ANe).* (R_Pr.ANe'*(W2_N .* (1 - J.Ne))) - J.ANe.*(R_An.ANe'*(W2_N .* J.Ne)));
        
        % diagonal terms
        FJAC = diag(1 + [Chi.Ne;Chi.ANe]) ...
            + diag([R_in_NES.Ne'*(W2_N .* J.Ne) + R_out_NES.Ne'*(W2_N .* (1-J.Ne)); ...
                  R_in_NES.ANe'*(W2_N .* J.ANe) + R_out_NES.ANe'*(W2_N .* (1-J.ANe))]) ...
            + diag([R_Pr.Ne'*(W2_N .* (1 - J.ANe)) + R_An.Ne'*(W2_N .* J.ANe); ...
                  R_Pr.ANe'*(W2_N .* (1 - J.Ne)) + R_An.ANe'*(W2_N .* J.Ne)]);
        % NES terms
        FJAC(iJNe,iJNe) = FJAC(iJNe,iJNe) - (R_in_NES.Ne').*((1 - J.Ne)*(W2_N)') ...
            - (R_out_NES.Ne').*(J.Ne*(W2_N)');
        FJAC(iJANe,iJANe) = FJAC(iJANe,iJANe) - (R_in_NES.ANe').*((1 - J.ANe)*(W2_N)') ...
            - (R_out_NES.ANe').*(J.ANe*(W2_N)');
        
        % Pair terms
        FJAC(iJNe,iJANe) = FJAC(iJNe,iJANe) + (R_Pr.Ne').*((1 - J.Ne)*(W2_N)') ...
            + (R_An.Ne').*(J.Ne*(W2_N)');
        FJAC(iJANe,iJNe) = FJAC(iJANe,iJNe) + (R_Pr.ANe').*((1 - J.ANe)*(W2_N)') ...
            + (R_An.ANe').*(J.ANe*(W2_N)');
        
        dJ = FJAC\FVEC;
        
        Jnew.Ne = J.Ne - dJ(iJNe);
        Jnew.ANe = J.ANe - dJ(iJANe);
        
    
        if ((norm(Jnew.Ne - J.Ne) <= Rtol * norm(J_in.Ne)) && (norm(Jnew.ANe - J.ANe) <= Rtol * norm(J_in.ANe)))
            CONVERGED_IN = true;
        end

        J = Jnew;

%         if (~CONVERGED_IN)
% 
%             eta_NES.Ne = R_in_NES.Ne' * (W2_N .* J.Ne);
%             Chi_NES.Ne = eta_NES.Ne + R_out_NES.Ne' * (W2_N .* (1 - J.Ne));
% 
%             eta_NES.ANe = R_in_NES.ANe' * (W2_N .* J.ANe);
%             Chi_NES.ANe = eta_NES.ANe + R_out_NES.ANe' * (W2_N .* (1 - J.ANe));
% 
%             eta_TP.Ne = R_Pr.Ne' * (W2_N .* (1 - J.ANe));
%             Chi_TP.Ne = eta_TP.Ne + R_An.Ne' * (W2_N .* J.ANe);
% 
%             eta_TP.ANe = R_Pr.ANe' * (W2_N .* (1 - J.Ne));
%             Chi_TP.ANe = eta_TP.ANe + R_An.ANe' * (W2_N .* J.Ne);
%         end
    end
    if(k_in >= maxIter)
        disp("Inner loop failed to converge within maxIter.");
    end

    iter_in = iter_in + k_in;
    
    % Update (U,J)
    Uvec(iY,mk) = 1 + C(iY) - Theta2_N' * (J.Ne - J.ANe)/ s_Y;
    Uvec(iE,mk) = 1 + C(iE) - Theta3_N' * (J.Ne + J.ANe) / s_E;
    
    Fvec(:,mk) = Uvec(:,mk) - U;
    
    if (mk == 1)
        Unew = Uvec(:,1);
    else
        alpha(1:mk-1) = (Fvec(:,1:mk-1) - Fvec(:,mk)*ones(1,mk-1))\(-Fvec(:,mk));
        alpha(mk) = 1 - sum(alpha(1:mk-1));
        Unew = Uvec(:,1:mk) * alpha(1:mk);
    end
    
    % check convergence
    if ((norm(Unew(1) - U(1)) <= Rtol * norm(U0(1))) && (norm(Unew(2) - U(2)) <= Rtol * norm(U0(2))))
        CONVERGED = true;
    end
    
    if (mk == m)
        Uvec = circshift(Uvec,-1,2);
        Fvec = circshift(Fvec,-1,2);
    end
    
    % update U, Y, E, T
    U = Unew;
    
    Y = U(iY) * Y0;
    E = U(iE) * E0 / Erg2MeV; % change unit    
    T = ComputeTemperatureFromSpecificInternalEnergy_TABLE(D, E, Y);
    
    % compute chemical potential and derivatives
    [Mnu, ~, ~] = ComputeNeutrinoChemicalPotentials(D, T, Y);
    
    % equilibrium distribution
    [J0.Ne, J0.ANe]= FermiDirac( g_E_N, Mnu, BoltzmannConstant * T );
    
    if (~CONVERGED)
        % FOR IMPLICIT R_NES, UPDATE (interpolate) R_NES here:
        [R_in_NES.Ne, R_out_NES.Ne, R_in_NES.ANe, R_out_NES.ANe] = ComputeNesScatteringOpacityOnEGrid_TABLE(g_E_N, D, T, Y);

        % scale R_NES by c * dt
        R_in_NES.Ne = R_in_NES.Ne * c * dt;
        R_out_NES.Ne = R_out_NES.Ne * c * dt;
        
        R_in_NES.ANe = R_in_NES.ANe * c * dt;
        R_out_NES.ANe = R_out_NES.ANe * c * dt;
        
        % update Pair kernels (R_pr and R_An)
        [R_Pr.Ne, R_An.Ne, R_Pr.ANe, R_An.ANe] = ComputePairOpacityOnEGrid_TABLE(g_E_N, D, T, Y);
        
        % scale R_Pr and R_An by c * dt
        R_Pr.Ne = R_Pr.Ne * c * dt;
        R_An.Ne = R_An.Ne * c * dt;
        
        R_Pr.ANe = R_Pr.ANe * c * dt;
        R_An.ANe = R_An.ANe * c * dt;
        
    end
end

if(k >= maxIter)
    disp("Failed to converge within maxIter.");
end


iter = k;
iter_in = iter_in / iter;


end












