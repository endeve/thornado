function [Jout, iter] = UpdateNeutrinoDistribution(Jin, J0, Chi)

maxIter = 10000;
Rtol = 1e-8;
% initial guess
Jout = Jin;

k = 0;
CONVERGED = false;

while((~CONVERGED)&&(k<=maxIter))
    k = k + 1;
    
    % Picard iteration
    Jnew = Jin + Chi.*(J0 - Jout);
    
    % check convergence
    if (norm(Jnew - Jout)< Rtol * norm(Jin))
        CONVERGED = true;
    end

    Jout = Jnew;
end
iter = k;

end

