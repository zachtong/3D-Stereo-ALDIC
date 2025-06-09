function cameraParams_matlab = cameraParamsFormatConvertFromMMC(CalibrationFilepath, CalibrationFile)
% This function is used to convert calibration format from MMC
% MMC is an open-source mutli-cameras calibration software:
% Refer to: https://github.com/JoeyYin-SEU/Multi_Camera_Calibration

cameraParams_MMC = load(fullfile(CalibrationFilepath, CalibrationFile),...
    'Camera_0001_Group_0001_IntrinsicMatrix','Camera_0001_Group_0001_RadialDistortion','Camera_0001_Group_0001_TangentialDistortion',...
    'Camera_0002_Group_0001_IntrinsicMatrix','Camera_0002_Group_0001_RadialDistortion','Camera_0002_Group_0001_TangentialDistortion',...
    'Camera_0002_Group_0001_R','Camera_0002_Group_0001_T');


cameraParams_matlab.cameraParamsLeft.K = cameraParams_MMC.Camera_0001_Group_0001_IntrinsicMatrix;
cameraParams_matlab.cameraParamsLeft.RadialDistortion = cameraParams_MMC.Camera_0001_Group_0001_RadialDistortion(1:3);
cameraParams_matlab.cameraParamsLeft.TangentialDistortion = cameraParams_MMC.Camera_0001_Group_0001_TangentialDistortion;
cameraParams_matlab.cameraParamsRight.K = cameraParams_MMC.Camera_0002_Group_0001_IntrinsicMatrix;
cameraParams_matlab.cameraParamsRight.RadialDistortion = cameraParams_MMC.Camera_0002_Group_0001_RadialDistortion(1:3);
cameraParams_matlab.cameraParamsRight.TangentialDistortion = cameraParams_MMC.Camera_0002_Group_0001_TangentialDistortion;
cameraParams_matlab.rotationMatrix = cameraParams_MMC.Camera_0002_Group_0001_R;
cameraParams_matlab.translationVector = cameraParams_MMC.Camera_0002_Group_0001_T;