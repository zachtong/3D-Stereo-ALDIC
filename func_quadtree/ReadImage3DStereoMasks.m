function [mask_file_name_0, mask_file_name_1, Img_0, Img_1, LoadImgMethod] = ReadImage3DStereoMasks(varargin)
%FUNCTION [file_name,Img,DICpara] = ReadImage3DStereo(varargin)
% MATLAB script: ReadImage3DStereo.m
% ----------------------------------------------
%   This script is to load 3D Stereo DIC images taken by two cameras. Images 
%   taken by Cameras 0 & 1 are named ended with '_0' and '_1', respectively.
%   Images can be loaded by:
%       i) selecting a folder which included all the DIC raw images, 
%       ii) inputing image file name prefix keywords
%       iii) manually select DIC raw images
%
%   INPUT: No inputs are needed
%
%   OUTPUT: file_name_0    Loaded DIC raw image file name taken by Camera 0
%           file_name_1    Loaded DIC raw image file name taken by Camera 1
%           Img_0          Loaded DIC images taken by Camera 0
%           Img_1          Loaded DIC images taken by Camera 1
%           DICpara      DIC parameters
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: aldicdvc@gmail.com, jin.yang@austin.utexas.edu
% Last time updated: 11/2023.
% ==============================================

%%
fprintf('Choose method to load image mask files:  \n')
fprintf('     0: Select image mask file folder;  \n')
fprintf('     1: Use prefix of image mask file names;  \n')
fprintf('     2: Manually select image mask files.  \n')
prompt = 'Input here: ';
LoadImgMethod = input(prompt);

switch LoadImgMethod 
    case 0
        % ==============================================
        imgfoldername = uigetdir(pwd,'Select image mask file folder');
        addpath([imgfoldername,'\']);
        %%% Camera 0 %%%
        img1 = dir(fullfile(imgfoldername,'*_0.jpg'));
        img2 = dir(fullfile(imgfoldername,'*_0.jpeg'));
        img3 = dir(fullfile(imgfoldername,'*_0.tif'));
        img4 = dir(fullfile(imgfoldername,'*_0.tiff'));
        img5 = dir(fullfile(imgfoldername,'*_0.bmp'));
        img6 = dir(fullfile(imgfoldername,'*_0.png'));
        img7 = dir(fullfile(imgfoldername,'*_0.jp2'));
        mask_file_name_0 = [img1;img2;img3;img4;img5;img6;img7];
        mask_file_name_0 = struct2cell(mask_file_name_0);

        %%% Camera 1 %%%
        img1 = dir(fullfile(imgfoldername,'*_1.jpg'));
        img2 = dir(fullfile(imgfoldername,'*_1.jpeg'));
        img3 = dir(fullfile(imgfoldername,'*_1.tif'));
        img4 = dir(fullfile(imgfoldername,'*_1.tiff'));
        img5 = dir(fullfile(imgfoldername,'*_1.bmp'));
        img6 = dir(fullfile(imgfoldername,'*_1.png'));
        img7 = dir(fullfile(imgfoldername,'*_1.jp2'));
        mask_file_name_1 = [img1;img2;img3;img4;img5;img6;img7];
        mask_file_name_1 = struct2cell(mask_file_name_1);

    case 1
        % ==============================================
        fprintf('What is prefix of DIC image mask files? E.g. img_0*.tif.   \n')
        prompt = 'Input here: ';
        mask_mask_file_name = input(prompt,'s');
        [~,imgname,imgext] = fileparts(mask_mask_file_name);
        mask_mask_file_name = dir([imgname,imgext]);
        mask_mask_file_name = struct2cell(mask_mask_file_name);

    otherwise
        % ==============================================
        disp('--- Please load first image mask file ---')
        mask_mask_file_name{1,1} = uigetfile('*.tif','Select reference Image (Deformed)');
        disp('--- Please load next image mask file ---')
        mask_mask_file_name{1,2} = uigetfile('*.tif','Select deformed Image (Reference)');
        prompt = 'Do you want to load more deformed image mask files? (0-Yes; 1-No)';
        DoYouWantToLoadMoreImages = input(prompt); imageNo = 2;
        while ( DoYouWantToLoadMoreImages == 0 )   
            imageNo = imageNo + 1;
            mask_mask_file_name{1,imageNo} = uigetfile('*.tif','Select next image mask file');
            prompt = 'Do you want to load more deformed image mask files? (0-Yes; 1-No)';
            DoYouWantToLoadMoreImages = input(prompt);
        end
end

% ==============================================
% Reduce 32/24 bits to 8 bits (Don't change original files)
numImages = size(mask_file_name_0,2);
for i = 1:numImages

    Img_0{i} = imread(mask_file_name_0{1,i});
    Img_1{i} = imread(mask_file_name_1{1,i});

    % Change color RGB images to grayscale images
    [~, ~, numberOfColorChannels_0] = size(Img_0{i});
    if numberOfColorChannels_0 == 4
        Img_0{i} = rgb2gray(Img_0{i}(:,:,1:3));  % ignore Alpha channel 
         
    elseif numberOfColorChannels_0 == 3
        Img_0{i} = rgb2gray(Img_0{i});
    
    end

    [~, ~, numberOfColorChannels_1] = size(Img_1{i});
    if numberOfColorChannels_1 == 4
         
        Img_1{i} = rgb2gray(Img_1{i}(:,:,1:3));  % ignore Alpha channel 
    elseif numberOfColorChannels_1 == 3
       
        Img_1{i} = rgb2gray(Img_1{i});
    end

    Img_0{i} = double(Img_0{i})';
    Img_1{i} = double(Img_1{i})';
    
end

% Change 32/24 bits to 8 bits (Change original files)
% ChangeImageIf32bits(file_name_0(1:2,:));
% ChangeImageIf32bits(file_name_1(1:2,:));

% ====== COMMENT ======
% Images physical world coordinates and image coordinates are different:
% --------------------
% --  This is image --
% |                  |
% y                  |
% |                  |
% |  --> x direction |
% |                  |
% --------------------
% after transforming,  MatLab matrix direction:
% --  This is matrix in Matlab --
% |                             |
% x                             |
% |                             |
% |  --> y direction            |
% |                             |
% --------------------------------
 
 
    
end

