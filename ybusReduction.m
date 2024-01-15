% Formation of Y Bus Matrix and Y Bus Reduction   Matlab Code
% Data given in an Excel file
% Author: Helal Uzzaman Hasib
% Date: 23 Oct 2023
% ===========================================================

% Clear the workspace and command window
clear
clc

% Data Loading from Excel File
problem = 'example_6.8';
fprintf('Given problem -->  "%s": \n', problem)
cd(problem);
A = xlsread('impedence_data');
cd ..
disp(A);


% Admittance matrix or Y-bus matrix formation
y = lfybus(A);y
n = length(y);

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
% %     y(d, : ) = [];  % remove specific row
% %     y(:, d) = [];   % remove specific column
%     disp('After removing specific bus:');
end

for i=1:length(RbusNo)
    d = RbusNo(i);
    y(d, : ) = [];  % remove specific row
    y(:, d) = [];   % remove specific column
end
disp(y);
