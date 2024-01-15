clear
clc

cd ('example_6.8p');
A = xlsread('impedence_data');
disp(A);
B = xlsread('bus_data');
disp(B);
cd ..
y = lfybusp(A);
N = size(y,1);

buses = B(:,1);
vmag = B(:,2);
delta = deg2rad(B(:,end));
v = vmag.*(cos(delta) + 1i*sin(delta));
pg = B(: , 3); qg = B(:, 4);
pl = B(: , 5); ql = B(:, 6); 
p = pg - pl;
q = qg - ql;
% bus selection 
pvbuses = buses( pg> 0 & buses ~= 1);
pqbuses = buses( pg <=0 & buses ~= 1); 
pvbuses = pvbuses';
pqbuses = pqbuses';


tol = 0.0001; maxIter = 100;

for iteration = 1:maxIter 
    % pqbuses 
    for i = pqbuses
        sum = 0;
        for n = 1: N
            if i~=n
                sum = sum + y(i,n)*v(n);
            end
        end
        s = p(i) - 1i*q(i);
        v(i) = 1/y(i,i) * (s/conj(v(i)) - sum);
    end
    % pvbuses
    for i = pvbuses
        I = 0;
        for n = 1: N
            I = I + y(i,n)*v(n);
        end
        q(i) = imag(v(i)*conj(I));
        sum = 0;
        for n = 1: N
            if i~=n
                sum = sum + y(i,n)*v(n);
            end
        end
        s = p(i) -1i*q(i);
        vnew = 1/y(i,i) * (s/conj(v(i)) - sum);
        a = angle(vnew);
        v(i) = abs(v(i)) * ( cos(a) + 1i* sin(a) );
    end
    vstore( iteration , :) = v;
    if iteration > 1
        diff = abs ( max(vstore(iteration,:)- vstore(iteration-1,:)));
        if diff< tol
            break
        end
    end
end

disp([ abs(vstore) rad2deg(angle(vstore)) ]);

