% 
% 
% Author: Helal Uzzaman Hasib
% Date: 17 Dec 2023
% ===========================================================

% Clear the workspace and command window
clear
clc

V.mag(1) = 10;  V.angle(1) = 90;
V.mag(2) =  5;  V.angle(2) = 30;
V.mag(3) = 10;  V.angle(3) = 210;

A.mag = 1;      A.angle = 120;
r = A.mag;
theta = deg2rad(A.angle);
a= r*(cos(theta) + 1i*sin(theta));
a

%  Matrix creating of original in Polar form
for k= 1: 3
    r = V.mag(k);
    theta = deg2rad( V.angle(k));
    v(k,1) = r*(cos(theta) + 1i*sin(theta));
end
v  % original matrix
% creating A matrix
AA = 1/3*[ 1  1    1;
      1  a  a^2;
      1 a^2   a;];
vaComponent = AA*v
AB = [ 1  1    1;
       1  a^2  a;
       1   a   a^2;];
  
vabcComponents= AB.*vaComponent
% original matrix ploting 
subplot(2,2,1)
compass(v)
title('ORIGINAL SYSTEM')
% gtext('SYMMETRICAL COMPONENTS')

subplot(2,2,2)
compass(vabcComponents(1,:));
title('ZERO SEQUENCE')
% gtext('SYMMETRICAL COMPONENTS')

subplot(2,2,3)
compass(vabcComponents(2,:));
title('POSITIVE SEQUENCE')
% gtext('SYMMETRICAL COMPONENTS')

subplot(2,2,4)
compass(vabcComponents(3,:));
title('NEGATIVE SEQUENCE')
% gtext('SYMMETRICAL COMPONENTS')
