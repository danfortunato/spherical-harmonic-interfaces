classdef sht_plan < handle
%SHT_PLAN   Spherical harmonic transform plan.
%   SHT_PLAN(DEG, N, BACKEND) creates a plan for a spherical harmonic
%   transform using all spherical harmonics up to degree DEG and an N x N
%   grid in latitude and longitude. Transforms between spherical harmonic
%   coefficients and values on the grid are done by the backend specified
%   by BACKEND. Currently, BACKEND may be:
%
%      'fasttransforms' - Mikael Slevinsky's FastTransforms library. [1]
%      'fmm3d' - The Flatiron Institute's FMM3D library. [2]
%      'shtns' - Nathanael Schaeffer's SHTns library. [3]
%
%   Each backend uses a different set of grid points on the sphere. The
%   exact grid points used by a plan are stored in PLAN.GRID.LAT (for
%   latitudinal grid points) and PLAN.GRID.LON (for longitudinal grid
%   points).
%
%   SHT_PLAN([LMAX MMAX], [NLAT NLON], BACKEND) uses spherical harmonics
%   with degree 0<=l<=LMAX and order |m|<=MMAX and an NLAT x NLON grid in
%   latitude and longitude.
%
%   [1] https://github.com/MikaelSlevinsky/FastTransforms
%   [2] https://github.com/flatironinstitute/FMM3D
%   [3] https://bitbucket.org/nschaeff/shtns

    properties

        lmax
        mmax
        nlat
        nlon
        grid
        method = 'shtns'

    end

    properties (Hidden)

        backend

    end

    methods

        function plan = sht_plan(varargin)

            narginchk(2, 3);

            if ( isscalar(varargin{1}) )
                [plan.lmax, plan.mmax] = deal(varargin{1});
            elseif ( isvector(varargin{1}) && length(varargin{1}) == 2 )
                plan.lmax = varargin{1}(1);
                plan.mmax = varargin{1}(2);
            else
                error('Bad input.');
            end

            if ( isscalar(varargin{2}) )
                [plan.nlat, plan.nlon] = deal(varargin{2});
            elseif ( isvector(varargin{2}) && length(varargin{2}) == 2 )
                plan.nlat = varargin{2}(1);
                plan.nlon = varargin{2}(2);
            else
                error('Bad input.');
            end

            if ( nargin == 3 )
                plan.method = varargin{3};
            end

            % Instantiate the backend
            switch lower(plan.method)
                case 'shtns'
                    plan.backend = backends.shtns.backend(plan);
                    plan.grid.lat = acos(util.legpts(plan.nlat));
                    plan.grid.lon = util.trigpts(plan.nlon, [0 2*pi]);
                case 'fmm3d'
                    plan.backend = backends.fmm3d.backend(plan);
                    plan.grid.lat = acos(util.legpts(plan.nlat));
                    plan.grid.lon = util.trigpts(plan.nlon, [0 2*pi]);
                case 'fasttransforms'
                    plan.backend = backends.fasttransforms.backend(plan);
                    plan.grid.lat = util.trigpts(plan.nlat, [0 pi]) + pi/(2*plan.nlat);
                    plan.grid.lon = util.trigpts(plan.nlon, [-pi pi]);
                otherwise
                    error('Unknown backend.');
            end

        end

        function V = coeffs2vals(plan, C)
            C = plan.backend.fromCanonicalCoeffs(C);
            V = plan.backend.coeffs2vals(C);
            V = plan.backend.toCanonicalVals(V);
        end

        function C = vals2coeffs(plan, V)
            V = plan.backend.fromCanonicalVals(V);
            C = plan.backend.vals2coeffs(V);
            C = plan.backend.toCanonicalCoeffs(C);
        end

        function C = toCanonicalCoeffs(plan, C)
            C = plan.backend.toCanonicalCoeffs(C);
        end

        function C = fromCanonicalCoeffs(plan, C)
            C = plan.backend.fromCanonicalCoeffs(C);
        end

        function V = toCanonicalVals(plan, V)
            V = plan.backend.toCanonicalVals(V);
        end

        function V = fromCanonicalVals(plan, V)
            V = plan.backend.fromCanonicalVals(V);
        end

    end

end
