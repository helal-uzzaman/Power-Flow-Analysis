% Formation of Y Bus Matrix and Y Bus Reduction   Matlab Code
% Data given in an Excel file
% Author: Helal Uzzaman Hasib
% Date: 23 Oct 2023
% ===========================================================

% Clear the workspace and command window
clear
clc
% cd('F:\1901017\power system\lab4')
A = xlsread('givenDataMam');
disp('Given data: ');
disp(A);
n = size(A,1);
for j=1:n
    z(A(j,2), A(j,3)) = A(j,4) + A(j,5)*i;
    z(A(j,3), A(j,2)) = A(j,4) + A(j,5)*i;
end
n = length(z);
for i=1:n
    for j=1:n
        if z(i,j) ==0
            z(i,j) = inf;
        end
    end
end
y = 1./z           
s = sum(y);
for i=1:n
    for j=1:n
        if i==j
            y(i,j)= s(i);
        else
            y(i,j) = -y(i,j);
        end
    end
end
disp('Y matrix: ');
disp(y);
% Y bus reduction 
n = length(y); % matrix length 
% n = 4;
% e = 2;
e = input('Total Number of bus has to be reduced: \n');
for r = 1:e
    d = input('Number of bus which has to be reduced:\n');
    n = length(y);
    for i= 1:n
        for j= 1: n
            y(i, j) = y(i,j) -(y(i,d)*y(d,j)/y(d,d));
        end
    end
    y(d, : ) = [];  % remove specific row
    y(:, d) = [];   % remove specific column
    disp(y);
end
% disp(Rbus);








