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
% cd('example_6.8');
A = xlsread('ex_impedence_data_3_bus');
disp('Given Impedance data from Excel file: ');      disp(A);


B = xlsread('ex_pv_data_3_bus');
disp('Load flow data from Excel file: ');             disp(B)
% cd ..

%===================== Y-Bus formation ======================

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

n = size(B, 1);
% formating complex power from given data
% load flow data extraction
PGCOL = 3; QGCOL = 4;
PLCOL = 5; QLCOL = 6;
% formating complex power from given data
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
    
    % calculate for load buses 
    for j=Load
        sum = 0;
        for k = 1: n
            if (j ~= k)
                sum = sum + Ybus(j, k)*v(k);
            end
        end
        % main formula implementation of Gauss-Seidel Method
        v(j) = 1/Ybus(j,j)* ((P(j)-Q(j))/ conj(v(j)) -sum);
    end 
    
    % calculation for generator (PV) buses

    
    for j=Gen
        % Q calculation 
        I = 0;
        for k=1: n
            I = I+ Ybus(j,k)*v(k);
        end
        Q(j) = 1i* imag(v(j)*conj(I));  % updating Q  working
        
        sum = 0;
        for k = 1: n
            if (j ~= k)
                sum = sum + Ybus(j, k)*v(k);
            end
        end
        % main formula implementation of Gauss-Seidel Method
  
        newV = 1/Ybus(j,j) * ( ((P(j)-Q(j))/ conj(v(j)) - sum));  %  working
        % Only imaginary part should be retained. |V| MUST HELD CONSTANT
        % SO REAL PART SHOULD BE CALCULATED.
        realPart = sqrt(abs(v(j))^2 - imag(newV)^2);    % working
        v(j) = realPart + i* imag(newV);     % working
        % using angle 
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
voltage
Vmag = abs(voltage);                Vangle = rad2deg(angle(voltage));
% % Output Decoration for Generalized
Iter = (1: iteration)';
for i=1: n
    eval(['V' num2str(i) '=Vmag(:,i);']);
    eval(['A' num2str(i) '=Vangle(:,i);']);
end
% T = table(Iter, V1, A1, V2, A2, V3, A3, V4, A4, V5 , A5);  
T = table(Iter, V1, A1, V2, A2, V3, A3);  
disp('Calculated Bus voltage|here V -> magnitude(volt) & A -> Angle(deg)|All in PU:'); disp(T);

