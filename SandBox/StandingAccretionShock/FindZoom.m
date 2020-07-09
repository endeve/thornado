% Spatial Grid:
eta = 0.5/8000.0;
N   = 512;

f = @(z) eta.*(z.^N-1.0)-(z-1.0);

format long
zoom = fzero( f, [ 1.0001, 1.8 ] )