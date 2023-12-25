% Formation of Y Bus Matrix and Y Bus Reduction   Matlab Code
% Data given in an Excel file
% Author: Helal Uzzaman Hasib
% Date: 23 Oct 2023
% ===========================================================

% Clear the workspace and command window
clear
clc
A = xlsread('givenDataYbusReduction');
disp('Given data in Excel: ');
disp(A);
n = size(A,1);
for j=1:n
    z(A(j,2), A(j,3)) = A(j,4) + A(j,5)*i;
    z(A(j,3), A(j,2)) = A(j,4) + A(j,5)*i;
end
n = length(z);
% Impedance matrix formation
for i=1:n
    for j=1:n
        if z(i,j) == 0
            z(i,j) = inf;
        end
    end
end
% Admittance matrix or Y-bus matrix formation
y = 1./z;          
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
e = input('Total Number of bus has to be reduced: \n');

for r = 1:e
    d = input('Number of bus which has to be reduced:\n');
    RbusNo(r) = d;
%     n = length(y);
    for i= 1:n
        for j= 1: n
            y(i, j) = y(i,j) -(y(i,d)*y(d,j)/y(d,d));  % used formula
        end
    end
%     y(d, : ) = [];  % remove specific row
%     y(:, d) = [];   % remove specific column
    disp('After removing specific bus:');
    disp(y);
end

for i=1:length(RbusNo)
    d = RbusNo(i);
    y(d, : ) = [];  % remove specific row
    y(:, d) = [];   % remove specific column
end
disp(y);
