% Symmetrical component visualisation
% Of a three phase system
% Author: Helal Uzzaman Hasib
% Date: 17 Dec 2023
% ===========================================================

% Clear the workspace and command window
clear
clc
% Input Matix
D = xlsread('systemData');
disp('Given three phase system data from Excel file: ');   disp(D);

V.mag(1) = D(1,1);  V.angle(1) = D(1,2);
V.mag(2) = D(2,1);  V.angle(2) = D(2,2);
V.mag(3) = D(3,1);  V.angle(3) = D(3,2);

% creating a(1<120) and a2(1<240)
A.mag = 1;      A.angle = 120;
r = A.mag;
theta = deg2rad(A.angle);
a= r*(cos(theta) + 1i*sin(theta));
a2 = a^2;

% creating original system in Complex form
v = zeros(3,1);
for k= 1: 3
    r = V.mag(k);
    theta = deg2rad( V.angle(k));
    v(k,1) = r*(cos(theta) + 1i*sin(theta));
end
disp('Unbalance system in complex form in a matrix: ')
disp(v') % original matrix
% creating zero, positive and negative sequence components of A phase
va0 = 1/3*(v(1) + v(2) + v(3));
va1 = 1/3*(v(1) + a*v(2) + a2*v(3));
va2 = 1/3*(v(1) + a2*v(2) + a*v(3));
% plotting the phasor of original and symmetrical Components of the System
vZeroSequence = [ va0     va0    va0]
vPosSequence  = [ va1  a2*va1  a*va1]
vNegSequence  = [ va2  a*va2  a2*va2]
 
% original matrix plotting 
subplot(2,2,1)
compass(v)
title('ORIGINAL SYSTEM')
% gtext('SYMMETRICAL COMPONENTS')

subplot(2,2,2)
compass(vZeroSequence);
title('ZERO SEQUENCE')
% gtext('SYMMETRICAL COMPONENTS')

subplot(2,2,3)
compass(vPosSequence);
title('POSITIVE SEQUENCE')
% gtext('SYMMETRICAL COMPONENTS')

subplot(2,2,4)
compass(vNegSequence);
title('NEGATIVE SEQUENCE')
% gtext('SYMMETRICAL COMPONENTS')
