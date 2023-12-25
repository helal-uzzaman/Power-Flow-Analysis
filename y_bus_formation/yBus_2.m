% Formation of Y Bus Matrix through Matlab Code
% Data given in an Excel file
% Author: Helal Uzzaman Hasib
% Date: 23 Oct 2023
% ===========================================================

% Clear the workspace and command window
clear
clc

% Read data from the Excel file named "givenData" into a matrix A
cd("/Users/helaluzzaman/Desktop/lab 3.2/power system/lab 4 ybusGen")
A = xlsread("givenData2");

% Display the given data
disp("Given data: ")
disp(A)

% Get the number of rows in the data
n = size(A,1);


% Construct impedance matrix from given data (excel data)
for i = 1:n
    % Extract data from the ith row of A and convert to complex impedance
    % Populate the impedance matrix Z with impedance values
    Z(A(i, 2), A(i, 3)) = A(i, 4) + 1i * A(i, 5);
    
    % Since the impedance matrix is symmetric, populate both Z(i,j) and Z(j,i)
    Z(A(i, 3), A(i, 2)) = A(i, 4) + 1i * A(i, 5);
end

% Get the number of rows in the impedance matrix
n = length(Z);

% Replace zero values in the impedance matrix with 'inf'
for i = 1:n
    for j = 1:n
        if Z(i, j) == 0
            Z(i, j) = inf;
        end
    end
end

% Calculate the admittance matrix Y by taking the reciprocal of Z
Y = 1./Z;

% Calculate the sum of admittances connected to each node
S = sum(Y);

% Populate the Y Bus matrix
for i = 1:n
    for j = 1:n
        if i == j
            % Setting the diagonal elements of Y Bus matrix
            Y(i, j) = S(i);
        else
            % Setting the off-diagonal elements of Y Bus matrix
            Y(i, j) = -Y(i, j);
        end
    end
end

% Display the Y Bus matrix
disp('Y-Bus matrix is : ');
disp(Y)
