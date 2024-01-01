% Load Flow Analysis using Gauss Seidel method
% Data given in two Excel file
% Author: Helal Uzzaman Hasib
% Date: 23 Oct 2023      Edited: 15 Dec 2023
% ===========================================================

% Clear the workspace and command window
clear
clc
maxIter = 100;      tol = 0.0001;
% Data Loading from Excel File
cd('example_6.7');
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



% n->  total bus number 
n = size(ybus,1);

for iter= 1: maxIter
    % calculation for pv buses
    for j = pvbuses
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
%         realPart = sqrt(abs(v(j))^2 - imag(vnew)^2);    % working
%         v(j) = realPart + i* imag(vnew);     % working
    end
    % calculation for pq buses 
    for j = pqbuses
        sum = 0;
        for k = 1: n
            if j ~=k
                sum = sum + ybus(j,k) * v(k);
            end           
        end
        sconj = (p(j) - 1i*q(j));
        v(j) = 1/ybus(j,j) *(sconj/conj(v(j)) - sum);
    end
    % storing the value
    voltage(iter, :) = v;
    
    % check iter > 1
    if iter > 1
        diff = abs ( max(voltage(iter,:)- voltage(iter-1,:)));
        if diff< tol
            break
        end
    end
end
voltagef = [abs(voltage) rad2deg(angle(voltage))]

Vmag = abs(voltage);                Vangle = rad2deg(angle(voltage));
% % Output Decoration for Generalized
Iter = (1: iter)';
for i=1: n
    eval(['V' num2str(i) '=Vmag(:,i);']);
    eval(['A' num2str(i) '=Vangle(:,i);']);
end
% T = table(Iter, V1, A1, V2, A2, V3, A3, V4, A4, V5 , A5);  
T = table(Iter, V1, A1, V2, A2, V3, A3);  
% disp('Calculated Bus voltage|here V -> magnitude(volt) & A -> Angle(deg)|All in PU:'); 
disp(T);





