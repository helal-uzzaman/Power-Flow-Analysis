function i = lfybusp(A)
N = size(A,1);

for i= 1:N
    R = A(i, 3);
    X = A(i, 4);
    row = A(i,1); col = A(i,2);
    z(row, col) = R + 1i*X;
    z(col, row) = R + 1i*X;
end
N = max(size(z));
z(z== 0) = inf;
y = 1./z; 
s = sum(y);
for i = 1: N
    
    for c = 1: N
        if (i == c)
            y(i,c) = s(i);
        else
            y(i,c) = -y(i,c);
        end
    end
end
i = y;

end