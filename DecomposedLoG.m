% Usage [xFil1, yFil1, xFil2, yFil2] = DecomposedLoG(hatRadius)
function [xFil1, yFil1, xFil2, yFil2] = DecomposedLoG(hatRadius)

% This function was developed by Lenny Campanello (ljcamp@umd.edu) in the
% lab of Wolfgang Losert (wlosert@umd.edu).  Please contact them for the
% correct citation to acknowledge if using this work.

% 2021/10/26 RML updated function name

cylRadiusFac = 3; 
kernelSize = cylRadiusFac*hatRadius;

x = -kernelSize:kernelSize;
y = -kernelSize:kernelSize;

sigR = hatRadius/sqrt(2); 

Gx = exp(-x.*x/2/sigR/sigR); Gx = Gx/sum(Gx(:)); 
Gy = exp(-y.*y/2/sigR/sigR); Gy = Gy/sum(Gy(:));  
Fx = (1/sigR/sigR)*(1 - x.*x/sigR/sigR); 
Fy = (1/sigR/sigR)*(1 - y.*y/sigR/sigR); 

xFil1 = Gx.*Fx; xFil1 = xFil1 - mean(xFil1); 
yFil1 = permute(Gy, [2, 1]); 

xFil2 = Gx; 
yFil2 = permute(Gy.*Fy, [2, 1]); yFil2 = yFil2 - mean(yFil2); 