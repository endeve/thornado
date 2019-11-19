function [MnuEff, TEff] = SearchEffectiveTandChemicalPotentials(J, E_N, iSpecies, M_matter, T_matter )

BoltzmannConstant =   8.6173e-11;

% y = kT
y0 = T_matter*BoltzmannConstant;

x0 = M_matter/y0;




x = linspace(x0-10,x0+10,100000);
% x = Mnu/kT;
if (iSpecies == 1)
    exponent = min( max( repmat(E_N,1,length(x)) - repmat(x,length(E_N),1), -100 ), 100 );
elseif (iSpecies == 2)
    exponent = min( max( repmat(E_N,1,length(x)) + repmat(x,length(E_N),1), -100 ), 100 );
else
    error('invalid species');
end
F = 1.0 ./ ( exp( exponent ) + 1 );
F2 = trapz(E_N, F.*E_N.^2);
F3 = trapz(E_N, F.*E_N.^3);

a = trapz(E_N, J.*E_N.^3);
b = trapz(E_N, J.*E_N.^2);

LHS = a./b.^(4/3);
RHS = F3./F2.^(4/3);


[val, idx] = min(abs(LHS - RHS));

eta = x(idx);
y2 = (b/F2(idx))^(1/3);
y3 = (a/F3(idx))^(1/4);

MnuEff = eta.*y2;
TEff = y2/BoltzmannConstant;


end
