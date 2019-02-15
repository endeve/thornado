function dF0dT = dFermiDiracdT(E_N, Mu, kT, dMudT, T )

    exponent = min( max( ( E_N - Mu ) / kT, -100 ), 100 );
    F0 = 1.0 ./ ( exp( exponent ) + 1 );
    dF0dT = ( F0 .* exp( exponent ) ) .* F0 .* ( dMudT + ( E_N - Mu ) / T ) / kT;

end
