% Newton raphson Load flow solution 

clear
clc
cd 'example_6.8'
A = xlsread('impedence_data');
B = xlsread('bus_data');
cd ..
y = lfybus(A);
B;
buses = B(:,1); 
theta = deg2rad(B(:,end));
v = B(:,2).*(cos(theta) + 1j*sin(theta)); v = v';
pg = B(:,3); 
p =  pg - B(:,5); p = p';
q = B(:,4) - B(:,6); q = q';

tol = 0.0001; maxIter = 100;
N = size(y,1);

% bus selection
slackbus = buses( buses==1);
pqbuses = buses( pg <= 0 & buses ~=1 )
pvbuses = buses(  pg > 0 & buses ~=1 )
buses = buses'; pqbuses = pqbuses' ; pvbuses = pvbuses';

% newton
Y = abs(y); theta = angle(y);
V = abs(v); delta = angle(v);


for iter = 1: maxIter
    pcal = zeros(1, N) ; qcal = zeros(1,N);
    for i= 2: N
        for n = 1: N
            pcal(i) = pcal(i) + V(i)*V(n)*Y(i,n)*cos(theta(i,n)- delta(i)+ delta(n));
            qcal(i) = qcal(i) - V(i)*V(n)*Y(i,n)*sin(theta(i,n)- delta(i)+ delta(n));
        end
    end
    delp = p - pcal; 
    delq = q - qcal;
    
    
    r = 1; % row initialisation 
    % J11 and J12 calculation 
    for i = [pqbuses pvbuses]
        % J11
        c = 1; % column initialisation 
        for j = [pqbuses pvbuses]
            if i~=j
                J(r,c) = -V(i)*V(j)*Y(i,j)*sin(theta(i,j)-delta(i)+delta(j));
            else
                J(r,c) = -qcal(i) - V(i)^2*Y(i,j)*sin(theta(i,j));
            end
            c = c + 1;
        end
        % J12
        for j = pqbuses
            if i~=j
                J(r,c) = V(i)*Y(i,j)*cos(theta(i,j)-delta(i)+delta(j));
            else
                J(r,c) = pcal(i)/V(i) + V(i)*Y(i,i)*cos(theta(i,i));
            end
            c = c + 1;
        end
        r = r+1;
    end
    % J21 and J22 calculation 
    for i = [pqbuses]
        % J21
        c = 1; % column re initialisation 
        for j = [pqbuses pvbuses]
            if i~=j
                J(r,c) = -V(i)*V(j)*Y(i,j)*cos(theta(i,j)-delta(i)+delta(j));
            else
                J(r,c) = pcal(i) - V(i)^2*Y(i,j)*cos(theta(i,j));
            end
            c = c + 1;
        end
        % J22
        for j = pqbuses
            if i~=j
                J(r,c) = -V(j)*Y(i,j)*sin(theta(i,j)-delta(i)+delta(j));
            else
                J(r,c) = qcal(i)/V(i) - V(i)*Y(i,j)*sin(theta(i,j));
            end
            c = c + 1;
        end
        r = r+1;
    end
    
    J
    
    F = [delp([pqbuses pvbuses]) delq(pqbuses)];
    
    
    F = F';
    delx = J\F;
    r = 1;  % initialisation of row for next assumption  calculation
    for i= [pqbuses pvbuses]
        delta(i) = delta(i) + delx(r);
        r = r+ 1;
    end
    for i= [pqbuses]
        V(i) = V(i) + delx(r);
        r = r+ 1;
    end
    rad2deg(delta)
    V
    % storing data
    voltage.mag(iter, :) = V; voltage.angle(iter, :)= delta;
    
    % convergence check
    if iter > 1
        if abs(max(delx)) < tol
            break
        end
    end
    
end

disp('result voltage magnitude and angle ');

iteration(:,1) = 1: iter;
t = table ( iteration, voltage.mag, rad2deg(voltage.angle) );
t.Properties.VariableNames = {'Iter', 'Voltage', 'Angle_degree'};
disp(t);


