function BicubicBspline = BicubicBsplineBuildFullTable(Img)
%%%%%%
%   Reference: Pan et al. https://doi.org/10.1016/j.taml.2016.04.003
%   Author: Zach
%%%%%%

BC = [-0.428571428571429	1.14285714285714	-1.14285714285714	0.428571428571429
1.01785714285714	-2.08928571428571	1.33928571428571	-0.267857142857143
-0.589285714285714	-0.0535714285714286	0.803571428571429	-0.160714285714286
0	1	0	0];

BicubicBspline.BicubicBsplineCoefTable = cell(size(Img));
BicubicBspline.BicubicBsplineCoefFlag = ones(size(Img));
BicubicBspline.FullOrNot = 1;

for i = 2:size(Img,1)-2
    for j = 2:size(Img,2)-2
        Q = Img(i-1:i+2,j-1:j+2);
        P = BC * Q * BC';
        A = rot90(P,2);
        BicubicBspline.BicubicBsplineCoefTable{i,j} = A;

    end
end


