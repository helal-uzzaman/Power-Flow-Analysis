clear
clc

cd ('example_6.8p');
A = xlsread('impedence_data');
disp(A);
cd ..

N = size(A,1);

for r= 1:N
    R = A(r, 3);
    X = A(r, 4);
    row = A(r,1); col = A(r,2);
    z(row, col) = R + 1i*X;
    z(col, row) = R + 1i*X;
end
N = max(size(z));
z(z== 0) = inf;
y = 1./z; 
s = sum(y);
for r = 1: N
    
    for c = 1: N
        if (r == c)
            y(r,c) = s(r);
        else
            y(r,c) = -y(r,c);
        end
    end
end
y

