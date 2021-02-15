function C = toPyramid(C)

if ( isvector(C) )
    n = sqrt(length(C))-1;
    if ( n ~= floor(n) )
        error('Coefficient vector must have length (n+1)^2.');
    end
else
    n = size(C,1)-1;
    m = (size(C,2)-1)/2;
    if ( n ~= m )
        error('Coefficient matrix must have size (n+1) x (2*n+1).');
    end
end

[i, j, idx] = util.pyramid(n);
if ( ~isvector(C) )
    C = C(idx);
end
C(abs(C) < eps) = 0;
A = zeros(n+1, 2*n+1);
A(idx) = C;
C = A;
%C = sparse(i, j, C, n+1, 2*n+1);

end
