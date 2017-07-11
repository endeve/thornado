close all
clear all

[ Time,~,~,~,~,~,E_Fk,~,~,~,~,~,~,~,~,~ ]...
  = ReadGlobalTally( 'LinearWaves1D_GlobalTally.dat' );

A = 1.0d-5;
tau = logspace(-1,10,11000)';

diff = zeros( size( tau, 1 ), 1 );
for i = 1 : size( tau, 1 )

   E_Fk_tau = 0.25 .* A.^2 .* cos(2 .* pi .* Time).^2 .* exp( - Time ./ tau(i) );
   diff(i) = norm( E_Fk_tau - E_Fk );
   
end

i_min_diff = find( diff == min( diff ) );
tau(i_min_diff)