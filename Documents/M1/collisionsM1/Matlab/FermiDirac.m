function F0 = FermiDirac(E_N, Mu, kT)

    exponent = min( max( ( E_N - Mu ) / kT, -100 ), 100 );
    F0 = 1.0 ./ ( exp( exponent ) + 1 );

end