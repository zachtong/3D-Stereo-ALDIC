function [DICpara] = setDICParas_STAQ(DICpara, file_name,Img,maskLeft,MaskRight,LoadImgMethod)
%SETDICPARAS_STAQ  Fill DIC subset parameters + derive ROI from mask.

% ROI from mask bounding box (accumulative) or full image (incremental)
if DICpara.DICIncOrNot == 0
    [rows, cols] = find(maskLeft{1});
    gridx(1) = min(rows);
    gridy(1) = min(cols);
    gridx(2) = max(rows);
    gridy(2) = max(cols);
    gridxy.gridx = gridx; gridxy.gridy = gridy;
else
    gridxy.gridx = [1,size(maskLeft{1},1)]; gridxy.gridy = [1,size(maskLeft{1},2)];
end



% Subset size
if isfield(DICpara, 'winsize') && ~isempty(DICpara.winsize)
    winsize = DICpara.winsize;
else
    fprintf('\n--- What is the subset size? ---\n');
    fprintf('Each subset has an area of [-winsize/2:winsize/2, -winsize/2:winsize/2]\n');
    winsize = input('Input an even number: ');
end

% Subset step
if isfield(DICpara, 'winstepsize') && ~isempty(DICpara.winstepsize)
    winstepsize = DICpara.winstepsize;
else
    fprintf('--- What is the subset step? ---\n');
    winstepsize = input('Input an integer to be a power of 2 (E.g., 16): ');
end

% ==============================================
% Subproblem 2 solver: only FEM is supported in the quadtree path
Subpb2FDOrFEM = 1;

% ==============================================
% Parallel cluster #
if isfield(DICpara, 'ClusterNo') && ~isempty(DICpara.ClusterNo)
    ClusterNo = DICpara.ClusterNo;
else
    ClusterNo = funParaInput('ClusterNo');
end


% ==============================================
if isfield(DICpara, 'winsizeMin') && ~isempty(DICpara.winsizeMin)
    winsizeMin = DICpara.winsizeMin;
else
    winsizeMin = funParaInput('winsizeMin');
end
DICpara.winsizeMin = winsizeMin;


% Inc-or-not, ImgSeqIncUnit, and NewFFTSearch are handled by
% setDICParas_IncOrNot.m (caller runs it first).

DICpara.winsize        = winsize;
DICpara.winstepsize    = winstepsize;
DICpara.gridxyROIRange = gridxy;
DICpara.LoadImgMethod  = LoadImgMethod;
DICpara.Subpb2FDOrFEM  = Subpb2FDOrFEM;
DICpara.ClusterNo      = ClusterNo;
DICpara.ImgSize        = size(Img{1});
DICpara.transformDisp  = [];

end

