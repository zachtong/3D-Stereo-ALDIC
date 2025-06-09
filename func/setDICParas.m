function [DICpara] = setDICParas(file_name,Img,LoadImgMethod) 
%% ============================================== 
% Decide DIC subset parameters 
numImages = size(file_name,2); 
 
% Choose ZOI 
 
fprintf('\n'); 
disp('--- Define ROI corner points at the top-left and the bottom-right ---') 
ImgTemp = imread(file_name{1}); 
if size(ImgTemp,3) > 3 
    ImgTemp = ImgTemp(:,:,1:3); % Delete alpha channel 
end 
imshow(ImgTemp);  
title('Click top-left and the bottom-right corner points','fontweight','normal','fontsize',16); 
 
gridx = zeros(1,2); gridy = zeros(1,2); 
[gridx(1), gridy(1)] = ginput(1); 
fprintf('Coordinates of top-left corner point are (%4.3f,%4.3f)\n',gridx(1), gridy(1)) 
 
[gridx(2), gridy(2)] = ginput(1); 
fprintf('Coordinates of bottom-right corner point are (%4.3f,%4.3f)\n',gridx(2), gridy(2)) 
 
gridxy.gridx = round(gridx); gridxy.gridy = round(gridy); 
 
% Choose subset size 
fprintf('\n'); 
fprintf('--- What is the subset size? --- \n'); 
fprintf('Each subset has an area of [-winsize/2:winsize/2, -winsize/2:winsize/2] \n'); 
prompt = 'Input an even number: '; 
winsize = input(prompt); 
 
% Choose subset size 
fprintf('--- What is the subset step? --- \n'); 
prompt = 'Input an integer: '; 
winstepsize = input(prompt); 
  
 
% ============================================== 
% Subproblem 2 solver: finite difference or finite element 
Subpb2FDOrFEM = 0; % By default initialize parameters 
Subpb2FDOrFEM = funParaInput('Subpb2FDOrFEM'); % Subproblem 2 using finite difference or fem? 
 
% ============================================== 
% Parallel cluster # 
ClusterNo = funParaInput('ClusterNo'); % Assign parpool cluster No 
  
% ============================================== 
% Deal with image sequence 
NewFFTSearch = 1; % By default initialize parameters 
if numImages > 2 
     
    % ============================================== 
    % DIC initial guess  
    NewFFTSearch = funParaInput('NewFFTSearch'); % Use last frame as init guess or not 
     
    % ============================================== 
    % Decide DIC as accumulative or incremental mode? 
    fprintf('--- Choose accumulative or incremental mode ---  \n') 
    fprintf('     0: Accumulative(By default);  \n') 
    fprintf('     1: Incremental;  \n') 
    prompt = 'Input here: '; 
    DICIncOrNot = input(prompt); 
     
    try 
        switch DICIncOrNot 
            case 0 
                ImgSeqIncUnit = numImages+1; 
                ImgSeqIncROIUpdateOrNot = 1; 
            case 1 
                % fprintf('Incremental mode: How many frames to update reference image once? \n'); 
                % prompt = 'Input here: '; 
                % ImgSeqIncUnit = input(prompt); 
                % fprintf('Update ROI at the same time of updating reference image? \n'); 
                % fprintf('    0: Do not update ROI; \n');  
                % fprintf('    1: Manually(Recommended); \n');  
                % fprintf('    2: Automatically; \n');  
                % prompt = 'Input here: '; 
                % ImgSeqIncROIUpdateOrNot = input(prompt); 
 
                % Default 
                ImgSeqIncUnit = 1; 
                ImgSeqIncROIUpdateOrNot = 0; 
 
            otherwise 
                ImgSeqIncUnit = numImages+1; 
                ImgSeqIncROIUpdateOrNot = 1; 
        end 
 
    catch 
        ImgSeqIncUnit = numImages+1;  
        ImgSeqIncROIUpdateOrNot = 1; 
    end 
 
     
     
% ================================     
else % Only two frames 
     
    ImgSeqIncUnit = numImages+1;  
    ImgSeqIncROIUpdateOrNot = 1; 
    DICIncOrNot = 0; % Revised by Zach  
end 
 
DICpara.winsize = winsize; 
DICpara.winstepsize = winstepsize; 
DICpara.gridxyROIRange = gridxy; 
DICpara.LoadImgMethod = LoadImgMethod; 
DICpara.DICIncOrNot = DICIncOrNot; % Revised by Zach 
DICpara.ImgSeqIncUnit = ImgSeqIncUnit; 
DICpara.ImgSeqIncROIUpdateOrNot = ImgSeqIncROIUpdateOrNot; 
DICpara.Subpb2FDOrFEM = Subpb2FDOrFEM; 
DICpara.NewFFTSearch = NewFFTSearch; 
DICpara.ClusterNo = ClusterNo; 
DICpara.ImgSize = size(Img{1}); 
 
end 
 
