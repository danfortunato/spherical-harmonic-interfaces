classdef backend < sht_backend

    properties
        z
        w
    end

    methods

        function backend = backend(plan)

            backend = backend@sht_backend(plan);

            if ( backend.mmax ~= backend.lmax )
                error('Discretization must be (P+1) x (2*P+1).');
            end

            [backend.z, backend.w] = util.legpts(backend.nlat);

        end

        function V = coeffs2vals(backend, C)
            nlat = backend.nlat;
            nlon = backend.nlon;
            z = backend.z;
            P = backend.lmax;
            P1 = P+1;
            P2 = 2*P+1;
            mex_id_ = 'shevalsphere(i dcomplex[xx], o dcomplex[xx], i int, i int, i int, i int, i double[x])';
[V] = backends.fmm3d.mex(mex_id_, C, P, P, nlat, nlon, z, P1, P2, nlon, nlat, nlat);
        end

        function C = vals2coeffs(backend, V)
            nlat = backend.nlat;
            nlon = backend.nlon;
            z = backend.z;
            w = backend.w;
            P = backend.lmax;
            P1 = P+1;
            P2 = 2*P+1;
            if ( isreal(V) )
                V = complex(V);
            end
            mex_id_ = 'projloc3d(i int, i int, i int, i int, i double[x], i double[x], i dcomplex[xx], o dcomplex[xx])';
[C] = backends.fmm3d.mex(mex_id_, P, P, nlat, nlon, z, w, V, nlat, nlat, nlon, nlat, P1, P2);
        end

        function [dlambda, dtheta] = coeffs2gradvals(backend, C)
            error('Not implemented.');
        end

        function C = toCanonicalCoeffs(backend, C)
            mmax = backend.mmax;
            C = [imag(C(:,1:mmax))*2*sqrt(2*pi) real(C(:,mmax+1))*sqrt(4*pi) real(C(:,mmax+2:end))*2*sqrt(2*pi)];
            C = util.fromPyramid(C);
        end

        function C = fromCanonicalCoeffs(backend, C)
            lmax = backend.lmax;
            mmax = backend.mmax;
            A = util.toPyramid(C);
            C = [zeros(lmax+1,mmax), A(:,mmax+1)/sqrt(4*pi), (A(:,mmax+2:end) - 1i*A(:,mmax:-1:1))/sqrt(2*pi)];
        end

        function V = toCanonicalVals(backend, V)
            V = real(V');
        end

        function V = fromCanonicalVals(backend, V)
            V = complex(V');
        end

    end

end
