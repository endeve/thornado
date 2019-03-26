function [Jout, iter] = UpdateNeutrinoDistribution_NES(Jin, J0, Chi, R_in_NES, R_out_NES)

global g_W2_N g_W3_N;

maxIter = 100;

% could possibly be looser
Rtol = 1e-8;

% initial guess
Jout = Jin;

k = 0;
CONVERGED = false;

W2_N = g_W2_N;
W3_N = g_W3_N;

while((~CONVERGED)&&(k<=maxIter))
    k = k + 1;
    
    % solve (66) and (75) in the tech report
    % Picard iteration
    % seems like the fixed-point has to be reformulated.
    
    % original formulation
%     NES_in = (1 - Jout) .* (R_in_NES'*diag(W2_N)*Jout);
%     NES_out = Jout .* (R_out_NES'*diag(W2_N)*(1 - Jout));
%     Jnew = (Jin + Chi.*J0 + NES_in - NES_out)./(1 + Chi);
    
    eta_NES = (R_in_NES'*diag(W2_N)*Jout);
    Chi_NES = (R_out_NES'*diag(W2_N)*(1 - Jout));
    
    
    % new formulation
    Jnew = (Jin + Chi.*J0 + eta_NES)./(1 + Chi + eta_NES + Chi_NES);
    
    % check detailed balance
    NES_in_E = (1 - Jout) .* (R_in_NES'*diag(W3_N)*Jout);
    NES_out_E = Jout .* (R_out_NES'*diag(W3_N)*(1 - Jout));
    norm(NES_in_E - NES_out_E)
    
    norm(Jnew - Jout)
    
    % check convergence
    if (norm(Jnew - Jout)< Rtol * norm(Jin))
        CONVERGED = true;
    end

    Jout = Jnew;
end
iter = k;

end

