% Spatial Grid:
eta = (1.8/256) / ( 1.8d0 );
N   = 128;

f = @(z) eta.*(z.^N-1.0)-(z-1.0);

format long
zoom = fzero( f, 1.2 )