classdef backend < sht_backend

    properties (Hidden)
        mwptr % A special property that MWrap uses as an opaque pointer to a C object
        mwptr_synthesis
        mwptr_analysis
    end

    methods

        function backend = backend(plan)

            backend = backend@sht_backend(plan);

            lmax = backend.lmax;
            mmax = backend.mmax;
            nlat = backend.nlat;
            nlon = backend.nlon;

            if ( mmax ~= lmax )
                error('Discretization must be (P+1) x (2*P+1).');
            end

            if ( nlat ~= lmax+1 )
                error('Number of coefficients and values must be the same.')
            end

            if ( nlon ~= 2*nlat-1 )
                error('Grid must be N x (2*N-1).');
            end

            n = backend.lmax+1;
            m = backend.nlon;

            % Create pointers to the C plans
            mex_id_ = 'o ft_harmonic_plan* = ft_plan_sph2fourier(i int)';
[mwptr] = backends.fasttransforms.mex(mex_id_, n);
            backend.mwptr = mwptr;

            mex_id_ = 'o ft_sphere_fftw_plan* = ft_plan_sph_synthesis(i int, i int)';
[mwptr_synthesis] = backends.fasttransforms.mex(mex_id_, n, m);
            backend.mwptr_synthesis = mwptr_synthesis;

            mex_id_ = 'o ft_sphere_fftw_plan* = ft_plan_sph_analysis(i int, i int)';
[mwptr_analysis] = backends.fasttransforms.mex(mex_id_, n, m);
            backend.mwptr_analysis = mwptr_analysis;

        end

        function delete(backend)
            if ~isempty(backend.mwptr)
                % Delete the C object
                mwptr_synthesis = backend.mwptr_synthesis;
                mwptr_analysis  = backend.mwptr_analysis;
                mex_id_ = 'ft_destroy_sphere_fftw_plan(i ft_sphere_fftw_plan*)';
backends.fasttransforms.mex(mex_id_, mwptr_synthesis);
                mex_id_ = 'ft_destroy_sphere_fftw_plan(i ft_sphere_fftw_plan*)';
backends.fasttransforms.mex(mex_id_, mwptr_analysis);
                mex_id_ = 'ft_destroy_harmonic_plan(i ft_harmonic_plan*)';
backends.fasttransforms.mex(mex_id_, backend);
                backend.mwptr           = '';
                backend.mwptr_synthesis = '';
                backend.mwptr_analysis  = '';
            end
        end

        % IMPORTANT! The polar grid is strange. The grids are:
        %    x = util.trigpts(nlon, [-pi pi]);
        %    y = util.trigpts(nlat, [0 pi]) + pi/(2*nlat);
        %    y = (0.5:nlat-0.5)'*pi/nlat;

        function V = coeffs2vals(backend, C)

            if ( ~isreal(C) )
                error('Coefficients must be real.')
            end

            n = backend.lmax+1;
            m = backend.nlon;
            ncoeffs = n*m;

            % Sphere to Fourier
            mex_id_ = 'ft_execute_sph2fourier(i ft_harmonic_plan*, io double[x], i int, i int)';
[C] = backends.fasttransforms.mex(mex_id_, backend, C, n, m, ncoeffs);

            % Fourier to grid
            mwptr_synthesis = backend.mwptr_synthesis;
            mex_id_ = 'ft_execute_sph_synthesis(i ft_sphere_fftw_plan*, io double[x], i int, i int)';
[C] = backends.fasttransforms.mex(mex_id_, mwptr_synthesis, C, n, m, ncoeffs);

            V = reshape(C, [n m]);

        end

        function C = vals2coeffs(backend, V)

            if ( ~isreal(V) )
                error('Values must be real.')
            end

            n = backend.lmax+1;
            m = backend.nlon;
            ncoeffs = n*m;

            % Grid to Fourier
            mwptr_analysis = backend.mwptr_analysis;
            mex_id_ = 'ft_execute_sph_analysis(i ft_sphere_fftw_plan*, io double[x], i int, i int)';
[V] = backends.fasttransforms.mex(mex_id_, mwptr_analysis, V, n, m, ncoeffs);

            % Fourier to sphere
            mex_id_ = 'ft_execute_fourier2sph(i ft_harmonic_plan*, io double[x], i int, i int)';
[V] = backends.fasttransforms.mex(mex_id_, backend, V, n, m, ncoeffs);

            C = reshape(V, [n m]);

        end

        function C = toCanonicalCoeffs(backend, C)
            mmax = backend.mmax;
            A = 0*C;
            A(:,mmax+1) = C(:,1);
            for m = 1:mmax
                A(m+1:end,mmax+1-m) = C(1:end-m,2*m);
                A(m+1:end,mmax+1+m) = C(1:end-m,2*m+1);
            end
            C = util.fromPyramid(A);
        end

        function C = fromCanonicalCoeffs(backend, C)
            mmax = backend.mmax;
            A = util.toPyramid(C);
            C = 0*A;
            C(:,1) = A(:,mmax+1);
            for m = 1:mmax
                C(1:end-m,2*m)   = A(m+1:end,mmax+1-m);
                C(1:end-m,2*m+1) = A(m+1:end,mmax+1+m);
            end
        end

        function V = toCanonicalVals(backend, V)
        end

        function V = fromCanonicalVals(backend, V)
        end

    end

end
