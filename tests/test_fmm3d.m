function pass = test_fmm3d()

tol = 1e-13;

% Test evaluation correctness on a random spherical harmonic expansion
rng(0)

% Size of coefficient matrix:
lmax = 20;   % Legendre
mmax = lmax; % Fourier

% Size of grid:
nlat = 200; % Latitude
nlon = 200; % Longitude

% Random spherical harmonic expansion:
c = randn((lmax+1)^2, 1);
c = sqrt(4*pi/nnz(c))*c; % Normalize so variance is 1

% FMM3D:
plan = sht_plan([lmax mmax], [nlat nlon], 'fmm3d');
V = plan.coeffs2vals(c);

% Direct:
A = util.sphvander(lmax, plan.grid.lat, plan.grid.lon);
U = A*c;

pass(1) = norm( U(:) - V(:), inf ) < tol;

% Now test that we get the same coefficients back:
d = plan.vals2coeffs(V);

pass(2) = norm( c - d, inf ) < tol;

end
