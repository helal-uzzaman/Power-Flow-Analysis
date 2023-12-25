% Load Flow Analysis using Newton Raphson method
% Data given in two Excel file
% Author: Helal Uzzaman Hasib
% Date: 26 Nov 2023
% ===========================================================

% Clear the workspace and command window
clear
clc
maxIter = 20;
tol = 0.001;

% Data Loading from Excel File
% cd('F:\1901017\power system\lab4')
A = xlsread('givenDataYbusExBook');
% A = xlsread('givenDataYbusLab');
disp('Given Impedance data from Excel file: ');      disp(A);

B = xlsread('loadFlowExBook');
% B = xlsread('loadFlowLab');
disp('Load flow data from Excel file: ');
disp(B)
% ===================== Y-Bus formation ======================
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

% =================== Newton Raphson Method ====================
% Ybus polar form conversion
Ymag = abs(Ybus);          % Ybus magnitude --> Ymag    
theta = angle(Ybus);      % Ybus angle     --> Yangle

n = size(B, 1);
% load flow data extraction
PGCOL = 3; QGCOL = 4;
PLCOL = 5; QLCOL = 6;
% formating complex power from given data
s = (B(:,PGCOL)+1i*B(:, QGCOL)) - (B(:,PLCOL)+1i* B(:,QLCOL));
p = real(s)';               % real power     --> P
q = imag(s)';               % reactive power --> q

Vangle = deg2rad(B(:,7));   % voltage angle --> Vangle
v = B(:,2).*(cos(Vangle) + 1i*sin(Vangle));   % taking initial voltage values
Vmag = abs(v);               % voltage magnitude --> Vmag

% iterate and find load and gen bus
n = size(B, 1);
l = 1; g = 1; Load = []; Gen = [];
for bus = 2: n
    if (B(bus, PGCOL) == 0 && B(bus, QGCOL)==0)   % checking if load bus 
        Load(l) = bus;
        l = l+1;
    elseif (B(bus, PGCOL) > 0 && B(bus, QGCOL)==0) % checking if gen bus
        Gen(g) = bus;
        g = g+1;
    else
        error('proper data not found for identify load and generator bus');
    end
end
fprintf('Load bus Number: '); disp(Load);
fprintf('Generator bus Number: '); disp(Gen);
fprintf('\n')
% iterate and find load and gen bus

%Calculation for initial power
n = length(p);
for k = 2: n
    pcal(k) = 0;
    qcal(k) = 0;
    for l = 1: n
        if k == l
            pcal(k) = pcal(k)+ (Vmag(k)^2)*Ymag(k,k)*cos(theta(k,k));
            qcal(k) = qcal(k) - (Vmag(k)^2)*Ymag(k,k)*sin(theta(k,k));
        else 

            pcal(k) = pcal(k) + Vmag(k)*Vmag(l)*Ymag(k,l)*cos(theta(k,l) + Vangle(l) - Vangle(k));
            qcal(k) = qcal(k) - Vmag(k)*Vmag(l)*Ymag(k,l)*sin(theta(k,l) + Vangle(l) - Vangle(k));
        end
    end
end

% finding delta power
delp = p(:)-pcal(:);
delq = q(:)-qcal(:);

% Finding Jacobian Matrix
totalBusN = length([Load Gen]);
loadBusN = length(Load);
J = zeros(2*loadBusN+length(Gen));

totalBus =[Load Gen];
Lbus = Load;

% H calculation H -> [totalBus X totalBus]
r = 1; % Row number intialization
for ii = totalBus % row
    c = 1; % column number intialization
    for jj = totalBus % col
        if ii == jj
            s = 0;       % summing all term then subtracting ii ==jj term
            for nn = 1: n
                s = s + Vmag(ii)*Vmag(nn)*Ymag(ii,nn)*sin(theta(ii,nn) + Vangle(nn) - Vangle(ii));
            end
            J(r, c) = s - Vmag(ii)*Vmag(jj)*Ymag(ii,jj)*sin(theta(ii,jj) + Vangle(ii) - Vangle(jj));
        else
            J(r, c) = -Vmag(ii)*Vmag(jj)*Ymag(ii,jj)*sin(theta(ii,jj) + Vangle(jj) - Vangle(ii));
        end
        c = c + 1;
    end
    r = r+1;
end
% L calculation L -> [totalBus X loadBus]
r = 1; % Row number intialization
for ii = totalBus
    c = totalBusN +1; % column number intialization
    for jj = Lbus
        if ii == jj
            s = 0;
            for nn = 1: n
                s = s + Vmag(nn)*Ymag(ii,nn)*cos(theta(ii,nn) + Vangle(nn) - Vangle(ii));
            end
            s = s - Vmag(jj)*Ymag(ii,jj)*cos(theta(ii,jj) + Vangle(jj) - Vangle(ii));
            J(r, c) = s + 2*Vmag(ii)*Ymag(ii, ii)*cos(theta(ii,jj));
        else
            J(r, c) = Vmag(ii)*Ymag(ii,jj)*cos(theta(ii,jj) + Vangle(jj) - Vangle(ii));
        end
        c = c+1;
    end
    r = r+1;
end
% M calculation M -> [loadBus X totalBus]
r = totalBusN + 1; % Row number intialization
for ii = Lbus % row
    c = 1; % column number intialization
    for jj = totalBus % col
        if ii == jj
            s = 0;
            for nn = 1: n
                s = s + Vmag(ii)*Vmag(nn)*Ymag(ii,nn)*cos(theta(ii,nn) + Vangle(nn) - Vangle(ii));
            end
            J(r, c) = s - Vmag(ii)*Vmag(jj)*Ymag(ii,jj)*cos(theta(ii,jj) + Vangle(ii) - Vangle(jj));
        else
            J(r, c) = -Vmag(ii)*Vmag(jj)*Ymag(ii,jj)*cos(theta(ii,jj) + Vangle(jj) - Vangle(ii));
        end
        c = c +  1;
    end
    r = r+1;
end
% N calculation N -> [loadBus X loadBus]
r = totalBusN + 1; % Row number intialization
for ii = Lbus % row
    c = totalBusN + 1; % column number intialization
    for jj = Lbus % col
        if ii == jj
            s = 0;
            for nn = 1: n
                s = s - Vmag(nn)*Ymag(ii,nn)*sin(theta(ii,nn) + Vangle(nn) - Vangle(ii));
            end
            s = s + Vmag(jj)*Ymag(ii,jj)*sin(theta(ii,jj) + Vangle(jj) - Vangle(ii));
            J(r, c) = s - 2*Vmag(ii)*Ymag(ii, ii)*sin(theta(ii,jj));
        else
            J(r, c) = -Vmag(ii)*Ymag(ii,jj)*sin(theta(ii,jj) + Vangle(jj) - Vangle(ii));
        end
        c = c + 1;
    end
    r = r+1;
end
fprintf('Jacobian Matrix: \n')
disp(J);

% delp and delq -> R matrix
r = 1;
for ii = totalBus
    C(r,1) = delp(ii);
    r = r+1;
end
for ii = Lbus
    C(r,1) = delq(ii);
    r = r+1;
end
fprintf('Value matrix: \n')
disp(C)

delR = J\C;
% Update data delta angle and voltage
r = 1;
for ii = totalBus
    Vangle(ii) = Vangle(ii) + delR(r);
    r = r+1;
end
for ii = Lbus
    Vmag(ii) = Vmag(ii) + delR(r);
    r = r+1;
end
Vangle' 
Vmag'



