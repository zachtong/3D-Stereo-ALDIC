function StereoInfo = StereoCameraCalibration()
% squareSize is the size of the squares on the calibration board (in meters).
% boardSize is the number of corners in each row and column of the calibration board.

% squareSize
% fprintf('What is the squareSize(mm)? \n')
% prompt = 'Input here: ';
% squareSize = input(prompt);

% Debug Zach
squareSize = 10;

% load images
[leftImagesInfo,rightImagesInfo,~,~,~] = ReadImage3DStereo;
leftImagesNames = fullfile(leftImagesInfo(2,:), leftImagesInfo(1,:));
rightImagesNames = fullfile(rightImagesInfo(2,:), rightImagesInfo(1,:));

% check if the images are 32 bits, if yes, then convert to 24 bits
ChangeImageIf32bits(leftImagesNames);
ChangeImageIf32bits(rightImagesNames);

calibrate_board_type = 'checkerboard';

switch calibrate_board_type
    case 'checkerboard' % Method 1: checkerboard
        % Detect corners of the calibration board
        [imagePoints,boardSize] = detectCheckerboardPoints(leftImagesNames,rightImagesNames);

        % debug
        worldPoints = generateCheckerboardPoints(boardSize,1);
        temp11 = worldPoints(:,1);
        worldPoints(:,1) = worldPoints(:,2);
        worldPoints(:,2) = temp11;

        for i = 1:size(imagePoints,3)
            temp_imagePoints = zeros(size(worldPoints,1),2);
            % 判断 - 改变顺序
            for j = 1:size(worldPoints,1)
                flag = [0 1;-1 0] * worldPoints(j,:)' + [ 0,sqrt(size(worldPoints,1))-1]';
                exchange_index = find(all(worldPoints == flag', 2) == 1);
                temp_imagePoints(exchange_index,:) = imagePoints(j,:,i,2);
            end
            % update index
            imagePoints(:,:,i,2) = temp_imagePoints;

            % a = imagePoints(:,:,i,1); b = imagePoints(:,:,i,2);
        end
    case 'circle_symmetric' % Method 2: circle_symmetric
        % Unfinished!!!

        % Input the dimension of symmetric circle grid. eg: 10*14
        dim_sym = [6,10];

        for i = 1:size(leftImagesNames, 2)
            % Left images
            temp_Image_L = imread(leftImagesNames{i});
            imagePoints(:,:,i,1) = detectCircleGridPoints(temp_Image_L,dim_sym,PatternType="symmetric");
            % Right images
            temp_Image_R = imread(rightImagesNames{i});
            imagePoints(:,:,i,2) = detectCircleGridPoints(temp_Image_R,[10,14],PatternType="symmetric");

        end
        worldPoints = generateCircleGridPoints(dim_sym,squareSize,PatternType="symmetric");
end

% % Check if there are enough valid images for calibration
% if size(imagesUsedLeft,1) < 2 || size(imagesUsedRight,1) < 2
%     error('Not enough valid images for calibration.');
% end

worldPoints = generateCheckerboardPoints(boardSize,squareSize);

% Calibrate camera paras
[cameraParams, ~, estimationErrors] = estimateCameraParameters(imagePoints, worldPoints);

% Output calibration results to the command window using fprintf
fprintf('Stereo Camera Calibration Results:\n');
fprintf('Left Camera Intrinsic Parameters:\n');
disp(cameraParams.CameraParameters1.K);
fprintf('Right Camera Intrinsic Parameters:\n');
disp(cameraParams.CameraParameters2.K);
fprintf('Rotation Matrix:\n');
disp(cameraParams.RotationOfCamera2);
fprintf('Translation Vector:\n');
disp(cameraParams.TranslationOfCamera2);
% fprintf('Average Reprojection Error: %.2f\n', estimationErrors);

% Store results in a structure
StereoInfo.cameraParams.cameraParamsLeft = cameraParams.CameraParameters1;
StereoInfo.cameraParams.cameraParamsRight = cameraParams.CameraParameters2;
StereoInfo.cameraParams.rotationMatrix = rotationMatrix;
StereoInfo.cameraParams.translationVector = translationVector;
StereoInfo.cameraParams.FundamentalMatrix = cameraParams.FundamentalMatrix;
StereoInfo.cameraParams.EssentialMatrix = cameraParams.EssentialMatrix;
StereoInfo.cameraParams.estimationErrors = estimationErrors;
% stereoCameras.reprojectionError = estimationErrors;
end
