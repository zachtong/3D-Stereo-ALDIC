function U = funSmoothDisp_Quadtree(U,DICmesh,DICpara)
%FUNSMOOTHDISP: U = funSmoothDisp(U,DICmesh,DICpara)
% Object: Smooth solved displacement field by gridfit
% ----------------------------------------------
%
%   INPUT: U                 Displacement vector: U = [Ux_node1, Uy_node1, Ux_node2, Uy_node2, ... , Ux_nodeN, Uy_nodeN]';
%          DICmesh           DIC mesh
%          DICpara           DIC parameters
%
%   OUTPUT: U                Smoothed displacement fields by gridfit
%
% ----------------------------------------------
% Reference
% [1] RegularizeNd. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/61436-regularizend
% [2] Gridfit. Matlab File Exchange open source. 
% https://www.mathworks.com/matlabcentral/fileexchange/8998-surface-fitting-using-gridfit
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 12/2020.
% ==============================================


%% Initialization
h = DICpara.winsizeMin;
coordinatesFEM = DICmesh.coordinatesFEM;
elementsFEM = DICmesh.elementsFEM;
winstepsize = DICpara.winstepsize;
% DispFilterSize = DICpara.DispFilterSize;
% DispFilterStd = DICpara.DispFilterStd;
DispFilterSize = 0;
DispFilterStd = 0;

FilterStd = DispFilterStd; FilterSizeInput = DispFilterSize; LevelNo = 1; 
 
         
%% prompt = 'Do you want to smooth displacement? (0-yes; 1-no)';
DoYouWantToSmoothOnceMore = 0; % DoYouWantToSmoothOnceMore = input(prompt);
if DoYouWantToSmoothOnceMore == 0  
    if isempty(FilterStd) == 1
        prompt = 'Choose filter standard deviation(0-default): ';
        FilterStd = input(prompt);
        if FilterStd == 0
            FilterStd = 0.5; 
        end
    else
        if FilterStd == 0
            FilterStd = 0.5;
        end
    end
    if isempty(FilterSizeInput) == 1
        prompt = 'Choose Gaussian filter size(0-default): ';
        FilterSizeInput = input(prompt);
        if FilterSizeInput == 0
            FilterSizeInput = 2*ceil(2*FilterStd)+1; 
        end
    else
        if FilterSizeInput == 0
            FilterSizeInput = 2*ceil(2*FilterStd)+1;
        end
    end
end

SmoothTimes = 1;
while (DoYouWantToSmoothOnceMore==0)
    Coordxnodes = [min(coordinatesFEM(:,1)):h:max(coordinatesFEM(:,1))]'; 
    Coordynodes = [min(coordinatesFEM(:,2)):h:max(coordinatesFEM(:,2))]';
    Iblur_Top10 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), U{1,1},Coordxnodes,Coordynodes,'regularizer','springs'); Iblur_Top10=Iblur_Top10';
    Iblur_Top20 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), U{1,2},Coordxnodes,Coordynodes,'regularizer','springs'); Iblur_Top20=Iblur_Top20';
    Iblur_Top30 = gridfit(coordinatesFEM(:,1), coordinatesFEM(:,2), U{1,3},Coordxnodes,Coordynodes,'regularizer','springs'); Iblur_Top30=Iblur_Top30';

    % Iblur_Top = nan(size(ULocal,1),1); Iblur_Top(2*CoordCrackTop-1) = ULocal(2*CoordCrackTop-1); Iblur_Top(2*CoordCrackTop) = ULocal(2*CoordCrackTop); 
    % Iblur_Top10 = reshape(Iblur_Top(1:2:end),M,N); Iblur_Top20 = reshape(Iblur_Top(2:2:end),M,N);
    % -------------------------------------------------------
    % Iblur_Top1 = reshape(imgaussfilt(Iblur_Top10,FilterStd,'FilterSize',FilterSizeInput,'FilterDomain','spatial'), M*N, 1);
    % Iblur_Top2 = reshape(imgaussfilt(Iblur_Top20,FilterStd,'FilterSize',FilterSizeInput,'FilterDomain','spatial'), M*N, 1);
    % -------------------------------------------------------
    imageFilter=fspecial('gaussian',FilterSizeInput,FilterStd);
    % Iblur_Top1 = reshape(nanconv(Iblur_Top10,imageFilter,'edge','nanout'), M*N, 1);
    % Iblur_Top2 = reshape(nanconv(Iblur_Top20,imageFilter,'edge','nanout'), M*N, 1);
    % ULocal(2*CoordCrackTop-1) = Iblur_Top1(CoordCrackTop); ULocal(2*CoordCrackTop) = Iblur_Top2(CoordCrackTop);
    Iblur_Top1 = nanconv(Iblur_Top10,imageFilter,'edge','nanout');
    Iblur_Top2 = nanconv(Iblur_Top20,imageFilter,'edge','nanout');
    Iblur_Top3 = nanconv(Iblur_Top30,imageFilter,'edge','nanout');
    for tempi = 1:size(coordinatesFEM,1)
        [row1,~] = find(Coordxnodes==coordinatesFEM(tempi,1));
        [row2,~] = find(Coordynodes==coordinatesFEM(tempi,2));

        U{1,1}(tempi) = Iblur_Top1(row1,row2);
        U{1,2}(tempi) = Iblur_Top2(row1,row2);
        U{1,3}(tempi) = Iblur_Top3(row1,row2);
    end
     
    % close all; Plotuv(ULocal,x,y); % Plotting u and v
    % close all; Plotdisp_show(ULocal,elementsFEM,coordinatesFEM);
     
    % prompt = 'Do you want to smooth displacement once more? (0-yes; 1-no)';
    
    SmoothTimes = SmoothTimes+1;
    if SmoothTimes > 2
        DoYouWantToSmoothOnceMore = 1;
    end
    % DoYouWantToSmoothOnceMore = input(prompt); 
    
end

