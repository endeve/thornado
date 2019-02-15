function [Jout, iter] = UpdateNeutrinoDistribution_exact(Jin, J0, Chi)

Jout = (Jin + Chi.*J0)./(1 + Chi);

iter = 1;

end

