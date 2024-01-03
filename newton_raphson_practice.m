% Load Flow Analysis using Newton Raphson method
% Data given in two Excel file
% Author: Helal Uzzaman Hasib
% Date: 15 Dec 2023
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
% same as gauss seidel 

% NR initialize data
n = length(ybus);
V = abs(v);     delta = angle(v);
Y = abs(ybus);  theta = angle(ybus);
pcal = zeros(1,n);
qcal = zeros(1,n);


for iter = 1: maxIter
    
    % calculation p and q
    for j = 2: n
        for k = 1: n
            pcal(j) = pcal(j) + V(j)*V(k)*Y(j,k)*cos(theta(j,k) - delta(j) + delta(k));
            qcal(j) = qcal(j) - V(j)*V(k)*Y(j,k)*sin(theta(j,k) - delta(j) + delta(k));
        end
    end
    % finding delta power
    delp = p(:)-pcal(:);
    delq = q(:)-qcal(:);
    
    % jacobian matrix formation 
    
    
    
    
    
    
    
end





