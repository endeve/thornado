function [Jout, iter] = UpdateNeutrinoDistribution(Jin, J0, Chi)

maxIter = 100;
Rtol = 1e-8;
% initial guess
Jout = Jin;

k = 0;
CONVERGED = false;

while((~CONVERGED)&&(k<=maxIter))
    k = k + 1;
    
    % Picard iteration
    % THIS IS NOT A CONTRACTION MAPPING (Chi_q>1 for some q)
    % DOES NOT CONVERGE
    Jnew = Jin + Chi.*(J0 - Jout);
    
    % check convergence
    if (norm(Jnew - Jout)< Rtol * norm(Jin))
        CONVERGED = true;
    end

    Jout = Jnew;
end
iter = k;

end

