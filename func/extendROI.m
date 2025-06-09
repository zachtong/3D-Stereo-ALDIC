function    [DICpara] = extendROI(DICpara,RD)
% For right image series, to extend the corresponding ROI to rectangle.

tempX = RD.ResultFEMesh_corr{1,1}(:, 1);
tempX(tempX > DICpara.ImgSize(1)-1 | tempX < 1) = nan;
tempY = RD.ResultFEMesh_corr{1,1}(:, 2);
tempY(tempY > DICpara.ImgSize(2)-1 | tempY < 1) = nan;

maxX = round(max(tempX)) + DICpara.winsize/2 + 1;
minX = round(min(tempX)) - DICpara.winsize/2 - 1;
maxY = round(max(tempY)) + DICpara.winsize/2 + 1;
minY = round(min(tempY)) - DICpara.winsize/2 - 1;

maxX = min([maxX, DICpara.ImgSize(1)-3]);
minX = max([minX, 4]);
maxY = min([maxY, DICpara.ImgSize(2)-3]);
minY = max([minY, 4]);

DICpara.gridxyROIRange.gridx = [minX, maxX];
DICpara.gridxyROIRange.gridy = [minY, maxY];

end
