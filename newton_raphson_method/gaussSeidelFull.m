% Load Flow Analysis using Gauss Seidel method
% Data given in two Excel file
% Author: Helal Uzzaman Hasib
% Date: 23 Oct 2023      Edited: 15 Dec 2023
% ===========================================================

% Clear the workspace and command window
clear
clc
maxIter = 100;      tol = 0.00001;
%===================== Y-Bus formation ======================
% cd('F:\1901017\power system\lab4')
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
B = xlsread('loadFlowLab');
disp('Load flow data from Excel file: ');             disp(B)
n = size(B, 1);
% formating complex power from given data
% load flow data extraction
PGCOL = 3; QGCOL = 4;
PLCOL = 5; QLCOL = 6;
% formating complex power from given data
s = (B(:,PGCOL)+1i*B(:, QGCOL)) - (B(:,PLCOL)+1i* B(:,QLCOL));
P = (B(:,PGCOL) - B(:,PLCOL));
Q = 1i *(B(:, QGCOL) -  B(:,QLCOL));

Vangle = deg2rad(B(:,7));   % voltage angle --> Vangle
v = B(:,2).*(cos(Vangle) + 1i*sin(Vangle));   % taking initial voltage values


% iterate and find load and gen bus
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

v = v';   



for iteration = 1: maxIter
    % Q calculation if any PV or Generator bus
    for ii= Gen
        sum = 0;
        for k=1: n
            sum = sum+ v(ii)*Ybus(ii,k)*v(k);
        end
        Q(ii) = i*imag(sum) ; % updating Q 
    end

    % calculate for gen buses
    
    for j=Gen
        sum = 0;
        for k = 1: n
            if (j ~= k)
                sum = sum + Ybus(j, k)*v(k);
            end
        end
        % main formula implementation of Gauss-Seidel Method
        newV = 1/Ybus(j,j)* (P(j)-Q(j)/ conj(v(j)) -sum);
        % Only imaginary part should be retained. |V| MUST HELD CONSTANT
        % SO REAL PART SHOULD BE CALCULATED.
        realPart = sqrt(abs(v(j))^2 - imag(newV)^2);
        v(j) = realPart + i* imag(newV);
        
    end 

    % calculate for load buses 
    for j=Load
        sum = 0;
        for k = 1: n
            if (j ~= k)
                sum = sum + Ybus(j, k)*v(k);
            end
        end
        % main formula implementation of Gauss-Seidel Method
        v(j) = 1/Ybus(j,j)* (P(j)-Q(j)/ conj(v(j)) -sum);
    end 
    voltage(iteration,:) = v;
    
    % check iter > 1
    if iteration > 1
        diff = abs ( max(voltage(iteration,:)- voltage(iteration-1,:)));
        if diff< tol
            fprintf('\nConverged At %d No. iteration. Tolerance: %f\n',iteration, tol);
            break
        end
    end
 
end
Vmag = abs(voltage);                Vangle = rad2deg(angle(voltage));
% Output Decoration for Generalized
Iter = (1: iteration)';
for i=1: n
    eval(['V' num2str(i) '=Vmag(:,i);']);
    eval(['A' num2str(i) '=Vangle(:,i);']);
end
T = table(Iter, V1,V2, V3, V4, V5, A1,A2,A3, A4, A5);   
disp('Calculated Bus voltage|here V -> magnitude(volt) & A -> Angle(deg)|All in PU:'); disp(T);

