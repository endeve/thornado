function [MnuEff, TEff] = ComputeEffectiveTandChemicalPotentials(J, E_N, iSpecies, M_matter, T_matter )

BoltzmannConstant =   8.6173e-11;

Tol = 1e-6;
if (iSpecies == 1)
    alpha = 1;
elseif (iSpecies == 2)
    alpha = -1;
else
    error('invalid species');
end

% y = kT
y0 = T_matter*BoltzmannConstant;

x0 = M_matter/y0;


Converge = 0;

% x = Mnu/kT;
x = x0;

% N = sum(g_W2_N .* J);
% E = sum(g_W3_N .* J);
N = trapz(E_N, J.*E_N.^2);
E = trapz(E_N, J.*E_N.^3);
c = N^(4/3)/E;

k = 0;
MaxIter = 10000;
while(~Converge && k < MaxIter)
    k = k + 1;
    if (iSpecies == 1)
        exponent = min( max( E_N - x, -100 ), 100 );
    elseif (iSpecies == 2)
        exponent = min( max( E_N + x, -100 ), 100 );
    else
        error('invalid species');
    end
    F = 1.0 ./ ( exp( exponent ) + 1 );
%     F2 = sum(g_W2_N .* F);
%     F3 = sum(g_W3_N .* F);
    F2 = trapz(E_N, F.*E_N.^2);
    F3 = trapz(E_N, F.*E_N.^3);
    G = F2^(4/3)/F3;
    x_new = x - alpha*(G - c);
    
    if (abs(x_new - x) < Tol)
        Converge = 1;
    end
    x = x_new;
end
if(k >= MaxIter)
    warning('max iteration');
end
if (iSpecies == 1)
    exponent = min( max( E_N - x, -100 ), 100 );
elseif (iSpecies == 2)
    exponent = min( max( E_N + x, -100 ), 100 );
else
    error('invalid species');
end
F = 1.0 ./ ( exp( exponent ) + 1 );
% F2 = sum(g_W2_N .* F);
F2 = trapz(E_N, F.*E_N.^2);

y = (N/F2)^(1/3);
MnuEff = x*y;
TEff = y/BoltzmannConstant;

end
