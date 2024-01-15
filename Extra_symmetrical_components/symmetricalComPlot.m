% symmetical components for 3 phase system
% Author: Helal Uzzaman Hasib
% Date: 15 Jan 2024
clear
clc

% system data change magnitude and angle as you wish
va = 10* (cosd(90) + 1i*sind(90));  % cosd() sind() takes degree as input
vb =  5* (cosd(30) + 1i*sind(30));
vc = 10* (cosd(210) + 1i*sind(210));
% system data change magnitude and angle as you wish

a = cosd(120) + 1i*sind(120);
a2 = a^2;

va0 = 1/3*(va + vb + vc);
va1 = 1/3*(va + a*vb + a2* vc);
va2 = 1/3*(va + a2*vb + a*vc);

vb0 = va0; vb1 = a2* va1; vb2 =  a*va2;
vc0 = va0; vc1 =  a* va1; vc2 = a2*va2;

subplot(2,2,1);
compass( [va vb vc]);
title('ORIGINAL SYSTEM')

subplot(2,2,2);
compass( [va0 vb0 vc0]);
title('ZERO SEQUENCE')

subplot(2,2,3);
compass( [va1 vb1 vc1]);
title('POSITIVE SEQUENCE')

subplot(2,2,4);
compass( [va2 vb2 vc2]);
title('NEGATIVE SEQUENCE')







