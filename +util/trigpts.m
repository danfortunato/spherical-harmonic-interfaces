function [x, w] = trigpts(n, dom)
%TRIGPTS   Equally-spaced nodes and weights.
%   [X, W] = TRIGPTS(N) returns N equally-spaced nodes X and weights W for
%   the periodic trapezoid rule on [-1,1).
%
%   [X, W] = TRIGPTS(N, [A B]) returns N equally-spaced nodes X and weights
%   W for the periodic trapezoid rule on [A,B).

if ( nargin == 0 )
    x = [];
    w = [];
    return
end

x = linspace(-1, 1, n+1).'; % Equally-spaced points
x = (x - x(end:-1:1)) / 2;  % Enforce symmetry
x(end) = [];                % Drop the right endpoint
w = 2/n * ones(1, n);       % Quadrature weights

if ( nargin > 1 )
    x = diff(dom)/2*(x+1) + dom(1);
    w = diff(dom)/2*w;
end

end
