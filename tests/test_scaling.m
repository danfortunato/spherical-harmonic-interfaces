%% Run a scaling benchmark on each backend
rng(0)

lmaxs = 2.^(4:12)-1;

t_shtns_plan               = zeros(numel(lmaxs), 1);
t_shtns_transform          = zeros(numel(lmaxs), 1);
t_fasttransforms_plan      = zeros(numel(lmaxs), 1);
t_fasttransforms_transform = zeros(numel(lmaxs), 1);
t_fmm3d                    = zeros(numel(lmaxs), 1);

k = 1;
for lmax = lmaxs

    % Fourier:
    mmax = lmax;

    % Size of grid:
    nlat = lmax+1;   % nlat = 2*(lmax+1);
    nlon = 2*nlat-1; % nlon = nlat;

    fprintf('%d x %d coefficients -> %d x %d values\n', lmax+1, 2*mmax+1, nlat, nlon);

    % Random spherical harmonic expansion
    C = randn((lmax+1)^2, 1);

    if ( lmax < 1024 )
        tic
        plan = sht_plan([lmax mmax], [nlat nlon], 'fmm3d');
        V = plan.coeffs2vals(C);
        t_fmm3d(k) = toc;
    end

    tic
    plan = sht_plan([lmax mmax], [nlat nlon], 'shtns');
    t_shtns_plan(k) = toc;
    tic
    V = plan.coeffs2vals(C);
    t_shtns_transform(k) = toc;

    tic
    plan = sht_plan([lmax mmax], [nlat nlon], 'fasttransforms');
    t_fasttransforms_plan(k) = toc;
    tic
    V = plan.coeffs2vals(C);
    t_fasttransforms_transform(k) = toc;

    k = k+1;

end

%% Plot
lw = 1;
fs = 16;

loglog(lmaxs, t_fmm3d,                    'ro-',  'LineWidth', lw), hold on
loglog(lmaxs, t_shtns_plan,               'bo--', 'LineWidth', lw)
loglog(lmaxs, t_shtns_transform,          'bo-',  'LineWidth', lw)
loglog(lmaxs, t_fasttransforms_plan,      'go--', 'LineWidth', lw)
loglog(lmaxs, t_fasttransforms_transform, 'go-',  'LineWidth', lw)
z = lmaxs(4:end);
loglog(z, 1e-8*z.^3, 'k--', 'LineWidth', lw)
hold off

title('$$(p+1) \times (2p+1)$$ spherical harmonics to values', 'Interpreter', 'Latex')
xlabel('$$p$$', 'interpreter', 'latex', 'FontSize', 1.2*fs)
ylabel('Time (s)', 'FontSize', fs)
legend('FMM3D', 'SHTns (plan)', 'SHTns (transform)', ...
    'FastTransforms (plan)', 'FastTransforms (transform)', 'Location', 'NorthWest')

axis tight
ylim([5e-4 20])
grid on
set(gca, 'FontSize', fs)
shg
