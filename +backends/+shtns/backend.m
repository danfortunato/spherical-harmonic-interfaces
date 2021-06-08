classdef backend < sht_backend

    properties
        mwptr % A special property that MWrap uses as an opaque pointer to a C object
        nspat
        nlm
    end

    methods

        function backend = backend(plan)

            backend = backend@sht_backend(plan);

            lmax = backend.lmax;
            mmax = backend.mmax;
            nlat = backend.nlat;
            nlon = backend.nlon;

            if ( mod(nlat, 2) ~= 0 )
                error('nlat must be even.')
            end

            % Create a pointer to the C plan
            mex_id_ = 'o shtns_info_t* = ishtns_init(i long, i long, i long, i long)';
[mwptr] = backends.shtns.mex(mex_id_, lmax, mmax, nlat, nlon);
            backend.mwptr = mwptr;

            % Store some quantities locally, for easier access
            mex_id_ = 'o long = get_nlm(i shtns_info_t*)';
[nlm] = backends.shtns.mex(mex_id_, backend);
            mex_id_ = 'o long = get_nspat(i shtns_info_t*)';
[nspat] = backends.shtns.mex(mex_id_, backend);
            backend.nlm = nlm;
            backend.nspat = nspat;

        end

        function delete(backend)
            if ~isempty(backend.mwptr)
                % Delete the C object
                mex_id_ = 'shtns_destroy(i shtns_info_t*)';
backends.shtns.mex(mex_id_, backend);
                backend.mwptr = '';
            end
        end

        function V = coeffs2vals(backend, C)
            nlm = backend.nlm;
            nspat = backend.nspat;
            mex_id_ = 'SH_to_spat(i shtns_info_t*, i dcomplex[x], o double[x])';
[V] = backends.shtns.mex(mex_id_, backend, C, nlm, nspat);
            V = reshape(V, [backend.nlat backend.nlon]);
        end

        function C = vals2coeffs(backend, V)
            nlm = backend.nlm;
            nspat = backend.nspat;
            mex_id_ = 'spat_to_SH(i shtns_info_t*, i double[x], o dcomplex[x])';
[C] = backends.shtns.mex(mex_id_, backend, V, nspat, nlm);
        end

        function [dlambda, dtheta] = coeffs2gradvals(backend, C)
            nlm = backend.nlm;
            nspat = backend.nspat;
            mex_id_ = 'SH_to_grad_spat(i shtns_info_t*, i dcomplex[x], o double[x], o double[x])';
[dtheta, dlambda] = backends.shtns.mex(mex_id_, backend, C, nlm, nspat, nspat);
            dlambda = reshape(dlambda, [backend.nlat backend.nlon]);
            dtheta  = reshape(dtheta,  [backend.nlat backend.nlon]);
        end

        function C = toCanonicalCoeffs(backend, C)
            lmax = backend.lmax;
            mmax = backend.mmax;
            mask = tril(true(lmax+1,mmax+1));
            A = zeros(lmax+1);
            A(mask) = real(C);
            B = zeros(lmax+1);
            B(mask) = -imag(C);
            C = [fliplr(B(:,2:end)) A];
            C = util.fromPyramid(C);
        end

        function C = fromCanonicalCoeffs(backend, C)
            lmax = backend.lmax;
            mmax = backend.mmax;
            A = util.toPyramid(C);
            C = [zeros(lmax+1,mmax), A(:,mmax+1), A(:,mmax+2:end) - 1i*A(:,mmax:-1:1)];
            mask = [false(lmax+1,mmax) tril(true(lmax+1,mmax+1))];
            C = C(mask);
        end

        function V = toCanonicalVals(backend, V)
            V = flipud(V);
        end

        function V = fromCanonicalVals(backend, V)
            V = flipud(V);
        end

    end

end
