function r = ybus(A)
n = size(A,1); % row number %  do not use length as it gives max(size(A))

for i = 1: n
    R = A(i,3);
    X = A(i,4);
    z(A(i,1), A(i,2)) = R + 1j*X;
    z(A(i,2), A(i,1)) = R + 1j*X; % for symmetry
end
z(z==0) = inf;
y = 1./z;

s = sum(y);
n = size(y,1);

for r = 1: n
    for c = 1: n
        if r == c
            y(r,c) = s(c);
        else
            y(r,c) = -y(r,c);
        end
    end
end
r = y;

end