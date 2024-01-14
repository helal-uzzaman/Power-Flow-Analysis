% Load Flow Analysis using Gauss Seidel method
% Data given in two Excel file
% Author: Helal Uzzaman Hasib
% Date: 23 Oct 2023      Edited: 31 Dec 2023
% ===========================================================

% Clear the workspace and command window
clear
clc
maxIter = 100;      tol = 0.0001;
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
v = v';
pg = B(:,3);  qg = B(:,4);
pl = B(:,5);  ql = B(:,6);
p = pg - pl;     q = qg - ql;
buses = B(:,1);

% bus selection
slackbus = buses(buses == 1);
pvbuses = buses(pg > 0 & buses ~=1);
pqbuses = buses(pg <= 0 & buses ~=1);
pvbuses = pvbuses'; pqbuses = pqbuses';   % important for -> for loop
fprintf('pq busses are : ');   disp(pqbuses);
fprintf('pv busses are : ');   disp(pvbuses);

n = size(ybus,1);     % n->  total bus number 
for iter= 1: maxIter
    for j = pvbuses            % calculation for pv buses
        % q calculation
        I = 0;
        for k = 1: n
            I = I + ybus(j,k)*v(k);
        end
        q(j) = imag(v(j)*conj(I));
        
        
        sum = 0; 
        for k = 1: n
            if j ~=k
                sum = sum + ybus(j,k) * v(k);
            end           
        end
        sconj = p(j) - 1i*q(j);
        vnew = 1/ybus(j,j) *(sconj/conj(v(j)) - sum);
        % only angle should be updated   
        theta = angle(vnew);
        v(j) = abs(v(j))*( cos(theta) + 1i* sin(theta)); 
    end
    for j = pqbuses          % calculation for pq buses 
        sum = 0;
        for k = 1: n
            if j ~=k
                sum = sum + ybus(j,k) * v(k);
            end           
        end
        sconj = (p(j) - 1i*q(j));
        v(j) = 1/ybus(j,j) *(sconj/conj(v(j)) - sum);
    end
    voltage(iter, :) = v;      % storing the value
    % Convergence check iter > 1
    if iter > 1
        diff = abs ( max(voltage(iter,:)- voltage(iter-1,:)));
        if diff< tol
            break
        end
    end
end

volt.mag = abs(voltage);
volt.angle = rad2deg(angle(voltage));

% % Output Decoration for Generalized
Iter = (1: iter)';
T = table(Iter, volt.mag, rad2deg(volt.angle));
T.Properties.VariableNames = {'Iter', 'Voltage', 'Angle_degree'};
disp(T);
disp('Voltage(Volt) are V1, V2, V3,... and Angle(Degree) are A1, A2, A3, ....')
disp('All in Per-Unit');
% Output Decoration for Generalized

