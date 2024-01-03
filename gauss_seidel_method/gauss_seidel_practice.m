% Load Flow Analysis using Gauss Seidel ONLY FOR PQ bus method
% Data given in two Excel file
% Author: Helal Uzzaman Hasib
% Date: 23 Oct 2023      Edited: 15 Dec 2023
% ===========================================================

% Clear the workspace and command window
clear
clc
maxIter = 10;      tol = 0.0001;
% Data Loading from Excel File
cd('example_6.8');
A = xlsread('impedence_data');
disp('Given Impedance data from Excel file: ');      disp(A);


B = xlsread('bus_data');
disp('Load flow data from Excel file: ');             disp(B)
cd ..
ybus = lfybus(A);
ybus

% bus data [ Bus Voltage Pg Qg Pl Ql angle ] extraction 

vmag = B(:,2); %voltage mag
vangle = deg2rad(B(: ,7)); % voltage in rad
v= vmag.*( cos(vangle) + 1i* sin(vangle)); 
v = v'
pg = B(:,3);  qg = B(:,4);
pl = B(:,5);  ql = B(:,6);
p = pg - pl;     q = qg - ql;

% for pq bus only system  6.7
n = size(ybus,1);

for iter= 1: maxIter
    
    for j = 2: n
        sum = 0;
        for k = 1: n
            if j ~=k
                sum = sum + ybus(j,k) * v(k);
            end           
        end
        sconj = p(j) - 1i*q(j);
        v(j) = 1/ybus(j,j) *(sconj/conj(v(j)) - sum);
        
        % storing the value
        voltage(iter, :) = v;
    end
end
voltagef = [abs(voltage) rad2deg(angle(voltage))]




