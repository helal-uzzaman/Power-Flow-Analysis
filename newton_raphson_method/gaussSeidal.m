% Load Flow Analysis using Gauss Seidel method
% Data given in two Excel file
% Author: Helal Uzzaman Hasib
% Date: 23 Oct 2023
% ===========================================================

% Clear the workspace and command window
clear
clc
%===================== Y-Bus formation ======================
% cd('F:\1901017\power system\lab4')
% A = xlsread('givenDataYbusLab');
A = xlsread('givenDataYbusLab');
disp('Given Impedance data from Excel file: ');      disp(A);
n = size(A,1);
for j=1:n
    z(A(j,2), A(j,3)) = A(j,4) + A(j,5)*1i;
    z(A(j,3), A(j,2)) = A(j,4) + A(j,5)*1i;
end
n = length(z);
for m=1:n
    for j=1:n
        if z(m,j) ==0
            z(m,j) = inf;
        end
    end
end
y = 1./z;           s = sum(y);
for m=1:n
    for j=1:n
        if m==j
            y(m,j)= s(m);
        else
            y(m,j)= -y(m,j);
        end
    end
end
disp('Y-Bus matrix: ');   disp(y);    Ybus  = y;
%=================== Gauss Seidel Method ====================
% B = xlsread('loadFlowLab');
B = xlsread('loadFlowLab');
disp('Load flow data from Excel file: ');
disp(B)
n = size(B, 1);
% formating complex power from given data
s = (B(:,3)+1i*B(:, 4)) - (B(:,5)+1i* B(:,6));
Vangle = deg2rad(B(:,7)); 
v = B(:,2).*(cos(Vangle) + 1i*sin(Vangle));   % taking initial voltage values
v = v';
s = s';

iteration= input('Enter the number of iteration of GaussSeidel method: ');
% iteration = 10;
for iter = 1: iteration
    for j=2:n
        sum = 0;
        for k = 1: n
            if (j ~= k)
                sum = sum + Ybus(j, k)*v(iter, k);
            end
        end
        % main formula implementation of Gauss-Seidel Method
        v(iter+1, j) = 1/Ybus(j,j)* (s(j)/ conj(v(iter,j)) -sum);
        v(iter+1, 1) = v(iter,1);   %Storing slack bus voltage unchanged
    end    
end
disp('Bus voltage in complex form (raw): ');
disp(v);       % display in Complex form
% Converting complex to polar form & display output 
Vmag = abs(v);                Vangle = rad2deg(angle(v));
Iteration = (0:iteration)'; 
V1 = Vmag(:,1);   V2 = Vmag(:,2);   V3 = Vmag(:,3);
A1 = Vangle(:,1); A2 = Vangle(:,2); A3 = Vangle(:,3);
disp('Calculated Bus voltage | here V -> magnitude(volt) & A -> Angle(deg) | All in PU: ');
t = table(Iteration, V1, V2, V3, A1, A2, A3 ); disp(t);


