function [results,BicubicBspline] = BicubicBsplineInterp(BicubicBspline,Img,X,Y)
%%%%%%
%   
%   Reference: Pan et al. https://doi.org/10.1016/j.taml.2016.04.003
%   Author: Zach
%%%%%%

results = zeros(size(X));
tempX = X(:);
tempY = Y(:);
deltaX = rem(tempY,1);
deltaY = rem(tempX,1);
VecX = [ones(size(tempX)) deltaX  deltaX.^2  deltaX.^3 ];
VecY = [ones(size(tempX)) deltaY  deltaY.^2  deltaY.^3 ];

floor_tempX = floor(tempX);
floor_tempY = floor(tempY);
all_index = sub2ind(size(BicubicBspline.BicubicBsplineCoefFlag),floor_tempY,floor_tempX);

% Interp
usedBicubicBsplineCoefTable = BicubicBspline.BicubicBsplineCoefTable(sub2ind(size(BicubicBspline.BicubicBsplineCoefFlag),floor_tempY,floor_tempX));

for i = 1:size(all_index,1)
        InterpCoefMatrix = usedBicubicBsplineCoefTable{i};
        [row,col] = ind2sub(size(results),i);
        results(row,col) = VecX(i,:) * InterpCoefMatrix * VecY(i,:)';
end

% for i = 1:size(all_index,1)
%         % InterpCoefMatrix = usedBicubicBsplineCoefTable{i};
%         [row,col] = ind2sub(size(results),i);
%         results(row,col) = VecX(i,:) * BicubicBsplineCoefTable(:,:,all_index(i)) * VecY(i,:)';
% end
