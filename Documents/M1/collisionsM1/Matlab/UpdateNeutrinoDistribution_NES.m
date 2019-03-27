function [Jout, iter] = UpdateNeutrinoDistribution_NES(Jin, J0, Chi, R_in_NES, R_out_NES)

global g_W2_N;

maxIter = 100;
Rtol = 1e-8;
% initial guess
Jout = Jin;

k = 0;
CONVERGED = false;

%W2_N = 4 * pi * g_W2_N;

while((~CONVERGED)&&(k<=maxIter))
    k = k + 1;
    
    % solve (66) and (75) in the tech report
    % Picard iteration
    NES_in = (4 * pi - Jout) .* (R_in_NES'*diag(g_W2_N)*Jout);
    NES_out = Jout .* (R_out_NES'*diag(g_W2_N)*(4 * pi - Jout));
    Jnew = (Jin + Chi.*J0 + NES_in - NES_out)./(1 + Chi);
    
    % check convergence
    if (norm(Jnew - Jout)< Rtol * norm(Jin))
        CONVERGED = true;
    end

    Jout = Jnew;
end
iter = k;

end

