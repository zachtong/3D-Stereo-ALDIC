function    [DICpara] = extendROI_inc(DICpara,RD,ImgSeqNum,camera0OrNot)
% For incremental mode, to extend the corresponding ROI to rectangle.
% This function can be customed for particular demands

if strcmp(camera0OrNot,'notCamera0') % For right Camera
    tempX = RD.ResultFEMesh_corr{ImgSeqNum-1,1}(:, 1);
    tempX(tempX > DICpara.ImgSize(1)-1 | tempX < 1) = nan;
    tempY = RD.ResultFEMesh_corr{ImgSeqNum-1,1}(:, 2);
    tempY(tempY > DICpara.ImgSize(2)-1 | tempY < 1) = nan;
else
    tempX = RD.CoorEachFrame{ImgSeqNum-1,1}(:, 1);
    tempX(tempX > DICpara.ImgSize(1)-1 | tempX < 1) = nan;
    tempY = RD.CoorEachFrame{ImgSeqNum-1,1}(:, 2);
    tempY(tempY > DICpara.ImgSize(2)-1 | tempY < 1) = nan;
end 

% maxX = round(max(tempX)) + DICpara.winsize/2 + 1;
% minX = round(min(tempX)) - DICpara.winsize/2 - 1;
% maxY = round(max(tempY)) + DICpara.winsize/2 + 1;
% minY = round(min(tempY)) - DICpara.winsize/2 - 1;

maxX = round(max(tempX)) + 1;
minX = round(min(tempX)) - 1;
maxY = round(max(tempY)) + 1;
minY = round(min(tempY)) - 1;

maxX = min([maxX, DICpara.ImgSize(1)-3]);
minX = max([minX, 4]);
maxY = min([maxY, DICpara.ImgSize(2)-3]);
minY = max([minY, 4]);

DICpara.gridxyROIRange.gridx = [minX, maxX];
DICpara.gridxyROIRange.gridy = [minY, maxY];

end
