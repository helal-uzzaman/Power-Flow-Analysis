% Formation of Y Bus Matrix through Matlab Code
% Author: Helal Uzzaman Hasib
% Date: 23 Dec 2023

clear;
clc
cd example;
A = xlsread('ex_impedence_data_4_bus');
disp('Given impedence data: ');
disp(A);
cd ..;
% z matix making
n = size(A, 1);
for r = 1: n
    R = A(r, 3);
    X = A(r, 4);
    z(A(r, 1), A(r, 2)) = R + 1j*X;
    z(A(r, 2), A(r, 1)) = R + 1j*X;
end
z
n = length(z);
for r = 1: n
    for c = 1: n
        if z(r,c) == 0
            z(r,c) = inf;
        end
    end
end

% Replace zeros with 'inf' 
% shortcut
% z(z == 0) = inf;

y = 1./z

s= sum(y)     % for diagonal calculation 

for r=1:n
    for c= 1:n
        if r==c
            y(r,r) = s(r);
        else
            y(r,c) = -y(r,c);
            
        end
    end
end
y

