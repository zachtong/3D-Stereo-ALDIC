function [full_path_0, full_path_1, masks_0, masks_1] = ReadImage3DStereo_STAQ_mask_inc(imageLeft,imageRight)
% ==============================================
% For inc mode, users have to import the exact number of masks as raw
% images. (Update every frame right now.)

  fprintf('\n How do you import masks?  \n')
    fprintf('     0: Select entire mask folder(*_0: left, *_1: right);  \n')
    fprintf('     1: Separately select the left and right mask folders;  \n')
    prompt = 'Input here: ';
    LoadMaskMethod = input(prompt);

    switch LoadMaskMethod
        case 0
            % ==============================================
            imgfoldername = uigetdir(pwd,'Select images folder');
            addpath([imgfoldername,filesep]);
            %%% Camera 0 %%%
            img1 = dir(fullfile(imgfoldername,'*_0.jpg'));
            img2 = dir(fullfile(imgfoldername,'*_0.jpeg'));
            img3 = dir(fullfile(imgfoldername,'*_0.tif'));
            img4 = dir(fullfile(imgfoldername,'*_0.tiff'));
            img5 = dir(fullfile(imgfoldername,'*_0.bmp'));
            img6 = dir(fullfile(imgfoldername,'*_0.png'));
            img7 = dir(fullfile(imgfoldername,'*_0.jp2'));
            full_path_0 = [img1;img2;img3;img4;img5;img6;img7];
            full_path_0 = struct2cell(full_path_0);

            %%% Camera 1 %%%
            img1 = dir(fullfile(imgfoldername,'*_1.jpg'));
            img2 = dir(fullfile(imgfoldername,'*_1.jpeg'));
            img3 = dir(fullfile(imgfoldername,'*_1.tif'));
            img4 = dir(fullfile(imgfoldername,'*_1.tiff'));
            img5 = dir(fullfile(imgfoldername,'*_1.bmp'));
            img6 = dir(fullfile(imgfoldername,'*_1.png'));
            img7 = dir(fullfile(imgfoldername,'*_1.jp2'));
            full_path_1 = [img1;img2;img3;img4;img5;img6;img7];
            full_path_1 = struct2cell(full_path_1);

        case 1
            % ==============================================
            leftfoldername = uigetdir(pwd,'Select LEFT camera folder');
            rightfoldername = uigetdir(pwd,'Select RIGHT camera folder');

            addpath([leftfoldername, filesep]);
            addpath([rightfoldername, filesep]);

            %%% Camera LEFT (full_path_0) %%%
            img1 = dir(fullfile(leftfoldername,'*.jpg'));
            img2 = dir(fullfile(leftfoldername,'*.jpeg'));
            img3 = dir(fullfile(leftfoldername,'*.tif'));
            img4 = dir(fullfile(leftfoldername,'*.tiff'));
            img5 = dir(fullfile(leftfoldername,'*.bmp'));
            img6 = dir(fullfile(leftfoldername,'*.png'));
            img7 = dir(fullfile(leftfoldername,'*.jp2'));
            full_path_0 = [img1; img2; img3; img4; img5; img6; img7];
            full_path_0 = struct2cell(full_path_0);

            %%% Camera RIGHT (full_path_1) %%%
            img1 = dir(fullfile(rightfoldername,'*.jpg'));
            img2 = dir(fullfile(rightfoldername,'*.jpeg'));
            img3 = dir(fullfile(rightfoldername,'*.tif'));
            img4 = dir(fullfile(rightfoldername,'*.tiff'));
            img5 = dir(fullfile(rightfoldername,'*.bmp'));
            img6 = dir(fullfile(rightfoldername,'*.png'));
            img7 = dir(fullfile(rightfoldername,'*.jp2'));
            full_path_1 = [img1; img2; img3; img4; img5; img6; img7];
            full_path_1 = struct2cell(full_path_1);
        otherwise
            warning('Unrecognized method.');
    end

    % For acc mode, no need to show below
    % --- Check how many images were found ---
    % num_left  = size(full_path_0, 2);
    % num_right = size(full_path_1, 2);
    %
    % fprintf('\nNumber of LEFT  masks (Camera 0): %d\n', num_left);
    % fprintf('Number of RIGHT masks (Camera 1): %d\n', num_right);
    %
    % if num_left ~= num_right
    %     warning('The number of images in LEFT and RIGHT is not the same!');
    % end

[masks_0,masks_1] = convertMasksTo8bit(full_path_0,full_path_1);
end


% ==============================================
% Reduce 32/24 bits to 8 bits (Don't change original files)
function [Img_0,Img_1] = convertMasksTo8bit(full_path_0,full_path_1)
numImages = size(full_path_0,2);
Img_0 = cell(numImages,1);
Img_1 = cell(numImages,1);

for i = 1:numImages

    Img_0{i} = imread(full_path_0{1,i});
    Img_1{i} = imread(full_path_1{1,i});

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


    Img_0{i} = logical((Img_0{i}))';
    Img_1{i} = logical((Img_1{i}))';
end
end
