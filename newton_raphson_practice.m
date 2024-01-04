% Load Flow Analysis using Newton Raphson method
% Data given in two Excel file
% Author: Helal Uzzaman Hasib
% Date: 15 Dec 2023
% ===========================================================

% Clear the workspace and command window
clear
clc
maxIter = 20;      tol = 0.0001;
% Data Loading from Excel File
cd('example_lab');
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
% v = v';
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

% F matrix initialization ( preallocation for speed only )
F = zeros(length([pqbuses pvbuses pqbuses]), 1);
% Jacobian matrix initialization 
J = zeros(length([pqbuses pvbuses pqbuses]));

for iter = 1: maxIter
    
    % calculation p and q ; have to initialize pcal and qcal 
    pcal = zeros(1,n);      qcal = zeros(1,n);
    for j = 2: n
        for k = 1: n
            pcal(j) = pcal(j) + V(j)*V(k)*Y(j,k)*cos(theta(j,k) - delta(j) + delta(k));
            qcal(j) = qcal(j) - V(j)*V(k)*Y(j,k)*sin(theta(j,k) - delta(j) + delta(k));
        end
    end
    % finding delta power
    delp = p(:)-pcal(:);
    delq = q(:)-qcal(:);
    % F matrix formation 
    F = [delp([pqbuses pvbuses]) ;delq(pqbuses)];
    
    % jacobian matrix formation 
    % J11 and J12 i, j for bus numbers and r, c for J matrix position
    r = 1; % Row number intialization
    for i = [pqbuses pvbuses]
        c = 1; % column number intialization
        for j = [pqbuses pvbuses]
            if i ~= j
                J(r, c) = -V(i)*V(j)*Y(i,j)*sin(theta(i,j)-delta(i)+delta(j));
            else
                J(r, c) = -qcal(i) - V(i)^2*Y(i,i) * sin(theta(i,i));
            end
            c = c+1;
        end
        % J12 
        for j = pqbuses
            if i ~= j
                J(r, c) = V(i)*Y(i,j)*cos(theta(i,j)-delta(i)+delta(j));
            else
                J(r, c) = pcal(i)/ V(i) + V(i)*Y(i,i) *cos(theta(i,i));
            end
            c = c+1;
        end
        r = r+1;
    end
    % J21 and J22    i, j for bus numbers and r, c for J matrix position 
    for i = pqbuses
        c = 1; % column number re-initialization for J21 and J22
        for j = [pqbuses pvbuses] % J21
            if i ~= j
                J(r, c) = -V(i)*V(j)*Y(i,j)*cos(theta(i,j)-delta(i)+delta(j));
            else
                J(r, c) = pcal(i) - V(i)^2*Y(i,i) * cos(theta(i,i));
            end
            c = c + 1;
        end
        % J22
        co = length([pqbuses pvbuses]);   % column offset for J22
        ro = length([pqbuses pvbuses]);   % row offset for J22
        
        for j = pqbuses
            if i ~= j
                J(r, c) = -V(i)*Y(i,j)*sin(theta(i,j)-delta(i)+delta(j));
            else
                J(r, c) = qcal(i)/ V(i) - V(i)*Y(i,i) *sin(theta(i,i));
            end
            c = c + 1;
        end
        r = r+1;
    end
    delx = J\F;
    r = 1;
    for i = [pqbuses pvbuses]
        delta(i) = delta(i)+ delx(r);
        r = r+1;
    end
    for i = pqbuses
        V(i) = V(i)+ delx(r);
        r = r+1;
    end
    % storing data
    voltage.mag(iter, :) = V; voltage.angle(iter, :)= delta;
    
    % convergence check
    if iter > 1
        if abs(max(delx)) < tol
            fprintf('\nConverged At %d No. iteration. Tolerance: 0.00001\n',iter);
            break
        end
    end

end

disp([voltage.mag rad2deg(voltage.angle)]);




