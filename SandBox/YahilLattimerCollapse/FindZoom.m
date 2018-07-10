% Spatial Grid:
eta = 1.0 / ( 1.d4 );
N   = 128;

f = @(z) eta.*(z.^N-1.0)-(z-1.0);

format long
zoom = fzero( f, 1.2 )