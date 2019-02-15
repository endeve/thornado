function dF0dY = dFermiDiracdY(E_N, Mu, kT, dMudY)

    exponent = min( max( ( E_N - Mu ) / kT, -100 ), 100 );
    F0 = 1.0 ./ ( exp( exponent ) + 1 );
    dF0dY = ( F0 .* exp( exponent ) ) .* F0 * dMudY / kT;

end