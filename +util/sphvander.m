function V = sphvander(lmax, lat, lon)
%SPHVANDER   Vandermonde matrix for spherical harmonics.
%   V = SPHVANDER(LMAX, LAT, LON) returns the Vandermonde matrix V of
%   spherical harmonics up to degree LMAX evaluated on a tensor product
%   grid with nodes LAT in latitude and LON in longitude. The result is a
%   matrix of size LENGTH(LAT)*LENGTH(LON) x (LMAX+1)^2 which maps
%   spherical harmonic coefficients to values on the given grid. The
%   coefficients, which are parametrized by (L,M), are ordered by the
%   degree L. For example, for LMAX = 2, the ordering should be:
%
%      (0,0), (1,-1), (1,0), (1,1), (2,-2), (2,-1), (2,0), (2,1), (2,2).
%
%   The values are ordered so that nodes in the same longitude are grouped
%   together.

nlat = length(lat);
nlon = length(lon);
nvals = nlat*nlon;
ncoeffs = (lmax+1)^2;
V = zeros(nvals, ncoeffs);

coslat = cos(lat(:)).';
lon = lon(:).';

% Handle the zero degree term separately since it's simple.
V(:,1) = 1/sqrt(4*pi) * ones(nvals,1);

k = 2;
for l = 1:lmax
    % Normalization terms for the associated Legendre functions.
    m = (0:l).';
    a = (-1).^m./sqrt((1+double(m==0))*pi);
    
    % Compute the associated Legendre functions of cos(th) (co-latitude)
    % We will do one for the positive (including zero) order associated
    % Legendre functions and one for the negative.
    Gp = a .* legendre(l, coslat, 'norm');
    Gn = Gp(2:end,:,:);

    % Permute rows and columns to do the tensor product computation.
    Gp = permute(Gp, [2 3 1]);
    Gn = permute(Gn, [2 3 1]);
    
    % Multiply the associated Legendre polynomials by the correct Fourier
    % modes in the longitude variable and sum up the results.
    Vp = bsxfun(@mtimes, Gp, permute(cos((0:l)'*lon), [3 2 1]));
    Vn = bsxfun(@mtimes, Gn, permute(sin((1:l)'*lon), [3 2 1]));
    Vp = reshape(Vp, [nvals l+1]);
    Vn = reshape(Vn, [nvals l]);

    % Storage order is [-m ... -1 0 1 ... m].
    V(:,k:k+2*l) = [fliplr(Vn) Vp];

    k = k+2*l+1;
end

end
