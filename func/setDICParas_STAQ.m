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



% Choose subset size
fprintf('\n');
fprintf('--- What is the subset size? --- \n');
fprintf('Each subset has an area of [-winsize/2:winsize/2, -winsize/2:winsize/2] \n');
prompt = 'Input an even number: ';
winsize = input(prompt);

% Choose subset size
fprintf('--- What is the subset step? --- \n');
prompt = 'Input an integer to be a power of 2 (E.g., 16): ';
winstepsize = input(prompt);


% ==============================================
% Subproblem 2 solver: only FEM is supported in the quadtree path
Subpb2FDOrFEM = 1;

% ==============================================
% Parallel cluster #
ClusterNo = funParaInput('ClusterNo'); % Assign parpool cluster No


% ==============================================
winsizeMin = funParaInput('winsizeMin'); % Assign the finest element size in the quadtree mesh
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

