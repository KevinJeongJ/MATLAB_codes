function [dLamH] = ascending2d(dLamH,x_min, x_max)

% this function is to sort the 2D matrix dLamH based on the values
% in the 1st row. based on the x_min and x_max, variables 
% not falling within the range are eliminated 
% 

% create a dummy matrix of the same dimensions as the dLamH
dummy = dLamH;

dummy = dummy(dummy(:,1)>x_min & dummy(:,1)<x_max,:);
% elinating rows for first column values that does not fall within the
% range

dummy = sortrows(dummy);
%sorting by the first column values

disp(dummy)
clear dLamH;

dLamH = dummy;
clear dummy;
end