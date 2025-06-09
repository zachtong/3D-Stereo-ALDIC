function BicubicBspline = BicubiclBsplineInit(Img)
%%%%%%
%
%   Reference: Pan et al. https://doi.org/10.1016/j.taml.2016.04.003
%   Author: Zach
%%%%%%


BicubicBspline.BicubicBsplineCoefTable = cell(size(ImgDef));
BicubicBspline.BicubicBsplineCoefFlag = zeros(size(ImgDef));

