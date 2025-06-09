function [localBicubicBsplineCoefTable] = BicubicBsplineBuild(Img,Coordinates)
%%%%%%
%   Reference: Pan et al. https://doi.org/10.1016/j.taml.2016.04.003
%   Author: Zach
%%%%%%

%% Coordinates' first element is the index of row in image matrix

% B = [-1 3 -3 1; 3 -6 3 0; -3 0 3 0; 1 4 1 0]/6;
% C = [71 -19 5 -1; -19 95 -25 5; 5 -25 95 -19; -1 5 -19 71]/56;
% BC = B*C;

BC = [-0.428571428571429	1.14285714285714	-1.14285714285714	0.428571428571429
1.01785714285714	-2.08928571428571	1.33928571428571	-0.267857142857143
-0.589285714285714	-0.0535714285714286	0.803571428571429	-0.160714285714286
0	1	0	0];


%localBicubicBsplineCoefTable = zeros(4,4,size(Coordinates,1));
localBicubicBsplineCoefTable = cell(size(Coordinates,1),1);
for num = 1:size(Coordinates,1)
i = Coordinates(num,1);% COL
j = Coordinates(num,2);% ROW

% Q = [Img(i-1,j-1) Img(i-1,j) Img(i-1,j+1) Img(i-1,j+2);
%     Img(i,j-1) Img(i,j) Img(i,j+1) Img(i,j+2);
%     Img(i+1,j-1) Img(i+1,j) Img(i+1,j+1) Img(i+1,j+2);
%     Img(i+2,j-1) Img(i+2,j) Img(i+2,j+1) Img(i+2,j+2)];
Q = Img(i-1:i+2,j-1:j+2);

P = BC * Q * BC';
A = rot90(P,2);
localBicubicBsplineCoefTable{num} = A;
end




