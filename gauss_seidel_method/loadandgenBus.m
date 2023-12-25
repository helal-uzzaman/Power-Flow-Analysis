clear
clc

B = xlsread('loadFlowExBook');
B = xlsread('loadFlowLab');
disp('Load flow data from Excel file: ');
disp(B)
% iterate and find load and gen bus
PGCOL = 3; QGCOL = 4;
PLCOL = 5; QLCOL = 5;

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