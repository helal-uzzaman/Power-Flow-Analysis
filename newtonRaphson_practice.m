clear
clc


cd( 'example_6.8p')
A = xlsread('impedence_data');
B = xlsread('bus_data');
cd ..
disp(A); disp(B);
y = lfybus(A);
% bus data extraction
buses = B(:,1);
V = B(:,2);
delta = deg2rad(B(:,end));
pg = B(:,3); qg = B(:,4); pl = B(:,5); ql = B(:,6);
p = pg -pl; q = qg - ql;
% bus selection 
pvbuses = buses( pg> 0 & buses ~= 1); pqbuses = buses ( pg<=0 & buses ~= 1);

% ybus in polar
Y = abs(y); theta = angle(y);
N = length (y); % bus number

tol = 0.0001; maxIter = 100;

for iter = 1: maxIter
    % cal p and g
    pcal = zeros(N, 1); qcal = zeros(N,1) ;
    for i = 2: N
        for n = 1: N
            pcal(i) = pcal(i) + V(i)*V(n)*Y(i,n)* cos(theta(i,n) - delta(i) + delta(n));
            qcal(i) = qcal(i) - V(i)*V(n)*Y(i,n)* sin(theta(i,n) - delta(i) + delta(n));
        end
    end
    delp = p- pcal;
    delq = q- qcal;
    F = [delp([pqbuses pvbuses ]) ;delq(pqbuses)];
    
    % jacobina matrix
    r = 1;
    for i = [pqbuses pvbuses]
        % J11
        c = 1;
        for j = [pqbuses pvbuses]
            if i~=j
                J(r,c) = - V(i)*V(j)*Y(i,j)* sin(theta(i,j) - delta(i) + delta(j));
            else
                J(r,c) = - qcal(i) - V(i)^2* Y(i,j)* sin(theta(i,j));
            end
            c = c+ 1;
        end
        % J12
        for j = pqbuses
            if i~=j
                J(r,c) =  V(i)*Y(i,j)* cos(theta(i,j) - delta(i) + delta(j));
            else
                J(r,c) = pcal(i)/V(i) + V(i)* Y(i,j)* cos(theta(i,j));
            end
            c = c+ 1;
        end
        r = r+1;
    end
    for i = pqbuses
        % J21
        c = 1;
        for j = [pqbuses pvbuses]
            if i~=j
                J(r,c) = - V(i)*V(j)*Y(i,j)* cos(theta(i,j) - delta(i) + delta(j));
            else
                J(r,c) = pcal(i) - V(i)^2* Y(i,j)* cos(theta(i,j));
            end
            
            c = c+ 1;
        end
        % J22
        for j = pqbuses
            if i~=j
                J(r,c) =  -V(i)*Y(i,j)* sin(theta(i,j) - delta(i) + delta(j));
            else
                J(r,c) = qcal(i)/V(i) - V(i)* Y(i,j)* sin(theta(i,j));
            end
            c = c+ 1;
        end
        r = r+1;
    end
    J
    delx = J\F;
    
    % new update
    
    r = 1;
    for i = [pqbuses pvbuses]
        delta(i) = delta(i) + delx(r);
        r = r+ 1;
    end
    for i = [pqbuses]
        V(i) = V(i) + delx(r);
        r = r+ 1;
    end
    % store data
    voltage.mag(iter, : ) = V;
    voltage.angle(iter, : ) = rad2deg(delta);
    if iter > 1
        diff = max(delx)
        if(diff< tol)
            break
        end
    end
end

disp( [voltage.mag voltage.angle])














