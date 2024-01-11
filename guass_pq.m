% gauss

% y bus formation
clear
clc
cd 'example_6.7'
A = xlsread('impedence_data');
B = xlsread('bus_data');
cd ..
y = lfybus(A);
B;
buses = B(:,1); buses = buses';
theta = B(:,end);
v = B(:,2).*(cos(theta) + 1j*sin(theta)); v = v';
p = B(:,3) - B(:,5); p = p';
q = B(:,4) - B(:,6); q = q';

tol = 0.0001; maxIter = 100;
N = size(y,1);

for iter = 1: maxIter
    % all pq buses
    for i = 2: N
        sum = 0;
        for n = 1:N
            if i~= n
                sum = sum + y(i,n)*v(n);
            end
        end
        v(i) = 1/y(i,i)*( (p(i)-1j*q(i))/conj(v(i)) - sum); 
    end
    
    vdata(iter, :) = v;
    % convergence check
    if iter > 1
        diff = abs( max(vdata(iter, : ) - vdata(iter-1, :)));
        if  diff< tol
            break;
        end
    end
    
end
disp('result voltage magnitude and angle ');

iteration(:,1) = 1: iter;
t = table ( iteration, abs(vdata), rad2deg(angle(vdata)) );
disp(t);


