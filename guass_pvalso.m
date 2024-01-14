% gauss

% y bus formation
clear
clc
cd 'example_6.8p'
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

for iter = 1: maxIter
    %  pv buses
    for i = pvbuses
        % Q calculation
        I = 0;
        for n = 1: N
            I = I + y(i,n)*v(n);
        end
        q(i) = imag(v(i) * conj(I));
        
        sum = 0;
        for n = 1:N
            if i~= n
                sum = sum + y(i,n)*v(n);
            end
        end
        vnew = 1/y(i,i)*( (p(i)-1j*q(i))/conj(v(i)) - sum); 
        theta = angle(vnew);
        v(i) = abs(v(i)) * (cos(theta) + 1j*sin(theta));
    end
    %  pq buses
    for i = pqbuses
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
t.Properties.VariableNames = {'Iter', 'Voltage', 'Angle_degree'};
disp(t);


