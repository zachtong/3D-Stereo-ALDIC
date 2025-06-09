image = imread(['C:\Users\13014\OneDrive - The University of Texas at Austin\Documents\MATLABCodes\StereoDIC_Challenge_1\StereoSample2 - Simulated\SimulatedTranslate\Minibatch\' ...
    'R\16-mm Step 04_1.tif']);
cameraParams = cameraParameters('IntrinsicMatrix', StereoInfo.cameraParams.cameraParamsRight.K',...
    'RadialDistortion', StereoInfo.cameraParams.cameraParamsRight.RadialDistortion);
undistortedImage = undistortImage(image, cameraParams);

% figure;
% subplot(1, 3, 1);
% imshow(image);
% title('Original Image');
% 
% subplot(1, 3, 2);
% imshow(undistortedImage);
% title('Undistorted Image');
% 
% subplot(1, 3, 3);
% imshow(image-undistortedImage);
% title('difference');

imwrite(undistortedImage,'C:\Users\13014\OneDrive - The University of Texas at Austin\Documents\MATLABCodes\StereoDIC_Challenge_1\StereoSample2 - Simulated\SimulatedTranslate\Minibatch\R_undistorted\16-mm Step 04_1.tif');