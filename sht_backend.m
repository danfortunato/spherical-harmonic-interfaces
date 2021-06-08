classdef (Abstract) sht_backend

    properties

        lmax
        mmax
        nlat
        nlon

    end

    methods

        function backend = sht_backend(plan)

            backend.lmax = plan.lmax;
            backend.mmax = plan.mmax;
            backend.nlat = plan.nlat;
            backend.nlon = plan.nlon;

        end

    end

    methods (Abstract)

        V = coeffs2vals(backend, C);
        C = vals2coeffs(backend, V);
        [dlambda, dtheta] = coeffs2gradvals(backend, C);

        C = toCanonicalCoeffs(backend, C);
        C = fromCanonicalCoeffs(backend, C);
        V = toCanonicalVals(backend, V);
        V = fromCanonicalVals(backend, V);

    end

end
