function C = fromPyramid(C)

n = size(C,1)-1;
m = (size(C,2)-1)/2;
if ( n ~= m )
    error('Coefficient matrix must have size (n+1) x (2*n+1).');
end

[~, ~, idx] = util.pyramid(n);
C = C(idx);

end
