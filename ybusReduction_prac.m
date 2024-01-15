clear
clc
cd ('example_6.8p');
A = xlsread('impedence_data');
disp(A);
cd ..
y = lfybusp(A);
N = length(y);
e = 1;
for b = 1: e
    n = 3; % bus number
    rbuses(b) = n;
    
    for i = 1: N
        for j = 1: N
            y(i,j) = y(i,j) - ( y(i,n)* y(n,j) /y(n,n))
        end
    end
end
y

for r = rbuses
    y ( :, r) = [];
    y (r,: ) = [];
end
y
    
