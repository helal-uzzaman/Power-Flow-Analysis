% y bus formation 
clear
clc
cd 'example_6.1'
A = xlsread('impedence_data');
cd ..
A

n = size(A,1); % row number %  do not use length as it gives max(size(A))

for i = 1: n
    R = A(i,3);
    X = A(i,4);
    z(A(i,1), A(i,2) ) = R + 1j*X;
    z(A(i,2), A(i,1) ) = R + 1j*X;
end
z(z==0) = inf;
y = 1./z;

s = sum(y);
s
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
y