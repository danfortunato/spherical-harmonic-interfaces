function [x, w] = legpts(n, dom)
%LEGPTS   Gauss-Legendre nodes and weights.
%   [X, W] = LEGPTS(N) returns N Gauss-Legendre nodes X and weights W for
%   Gauss-Legendre quadrature on [-1,1].
%
%   [X, W] = LEGPTS(N, [A B]) returns N Gauss-Legendre nodes X and weights
%   W for Gauss-Legendre quadrature on [A,B].

if ( nargin == 0 )
    x = [];
    w = [];
    return
end

% The Golub-Welsch algorithm:
beta = .5./sqrt(1-(2*(1:n-1)).^(-2)); % Three-term recurrence
T = diag(beta, 1) + diag(beta, -1);   % Jacobi matrix
[V, D] = eig(T);                      % Eigenvalue decomposition
x = diag(D);                          % Legendre points
[x, i] = sort(x);                     % Sort
w = 2*V(1,i).^2;                      % Quadrature weights

% Enforce symmetry:
x = x(1:floor(n/2));
w = w(1:floor(n/2));
if ( mod(n, 2) == 1 )
    x = [ x ; 0          ; -x(end:-1:1) ];
    w = [ w , 2-sum(2*w) ,  w(end:-1:1) ];
else
    x = [ x ; -x(end:-1:1) ];
    w = [ w ,  w(end:-1:1) ];
end

if ( nargin > 1 )
    x = diff(dom)/2*(x+1) + dom(1);
    w = diff(dom)/2*w;
end

end
