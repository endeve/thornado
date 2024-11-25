function [F0_Ne, F0_ANe] = FermiDirac(E_N, Mu, kT)

    exponent = min( max( ( E_N - Mu ) / kT, -100 ), 100 );
    F0_Ne = 1.0 ./ ( exp( exponent ) + 1 );
    
    % assuming that the positron chemical potential is negative the
    % electron chemical potential
    exponent = min( max( ( E_N + Mu ) / kT, -100 ), 100 );
    F0_ANe = 1.0 ./ ( exp( exponent ) + 1 );

end