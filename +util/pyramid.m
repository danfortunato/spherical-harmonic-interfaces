function [i, j, idx] = pyramid(n, trans)
%PYRAMID   Pyramid indexing.
%   [I, J, IDX] = PYRAMID(N) returns the row and column indices (I,J)
%   corresponding to the "pyramid entries" of a (N+1) x (2*N+1) matrix. The
%   indices are given in row-major order. The vector IDX contains the
%   linear indices corresponding to (I,J).
%
%   [I, J, IDX] = PYRAMID(N, true) returns the row and column indices (I,J)
%   corresponding to the "pyramid entries" of a (N+1) x (N+1) matrix. The
%   indices are given in column-major order. The vector IDX contains the
%   linear indices corresponding to (I,J).
%
%   Pyramid indexing is useful for storing spherical harmonic coefficients.
%
%   Example:
%
%      A = [0 0 1 0 0;
%           0 2 3 4 0;
%           5 6 7 8 9];
%      [i, j, idx] = pyramid(2);
%      norm(A(idx) - (1:9).')
%
%      A = [1 0 0;
%           2 4 0;
%           3 5 6];
%      [i, j, idx] = pyramid(2, 'col');
%      norm(A(idx) - (1:6).')

if ( nargin < 2 )
    trans = false;
end

if ( ~trans )
    I = repmat(1:n+1, 2*n+1, 1);
    J = repmat((1:2*n+1)', 1, n+1);
    mask = (I+J > n+1) & (J-I < n+1);
    [j, i] = find(mask);
else
    I = repmat(1:2*n+1, n+1, 1);
    J = repmat((1:n+1)', 1, 2*n+1);
    mask = (I+J > n+1) & (I-J < n+1);
    [i, j] = find(mask);
end

idx = mysub2ind([n+1 2*n+1], i, j);

end

function idx = mysub2ind(siz, i, j)
    idx = i + (j-1)*siz(1);
end